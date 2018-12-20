require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/Bellmunt'

Workflow.require_workflow "HTS"
Workflow.require_workflow "Sequence"
Workflow.require_workflow "Sample"
Workflow.require_workflow "HLA"
Workflow.require_workflow "CFDTL"

module Bellmunt
  extend Workflow
  def self.reference
    'b37'
  end

  def self.organism
    'feb2014'
  end


  def self.path
    Rbbt.share.projects.Bellmunt
  end

  def self.interval_list
		Rbbt.share.projects.Bellmunt.data['Padded.b37.bed'].find
  end

  def self.pon
		Rbbt.share.projects.Bellmunt.data['pon.vcf'].find
  end

  def self.af_not_in_resource
    0.0000025
  end

  def self.germline_resource
    GATK.known_sites.b37["af-only-gnomad.vcf.gz"].produce.find
	end

  dep_task :reliable_mutations, Bellmunt, :genomic_mutations

  extension :bam
  dep HTS, :BAM_rescore, 
    :fastq1 => :placeholder, :fastq2 => :placeholder, :reference => :placeholder,
    :sample_name => :placeholder,
    :platform_unit => :placeholder,
    :read_group_name => :placeholder,
    :sequencing_center => "CNAG",
    :platform => 'Illuimna',
    :library_name => 'LN',
    :interval_list => Bellmunt.interval_list do |jobname,options|

    sample_name = jobname

    options[:sample_name] = sample_name.gsub('_', '.')

    options[:reference] = Bellmunt.reference
    sample_fastqs = []
    Bellmunt.path.glob("*.fastq.gz").each do |file|
      basename = File.basename(file)
      sample_fastqs << file if basename.split("_")[2] == sample_name
    end

    case sample_fastqs.length 
    when 2
      machine, lane, sample = File.basename(sample_fastqs.first).split("_")
      options[:read_group_name] = [machine, lane] * "."
      options[:platform_unit] = [machine, lane, sample] * "."
      options[:fastq1] = sample_fastqs.sort.first
      options[:fastq2] = sample_fastqs.sort.last
      {:inputs => options, :jobname => jobname}
    when 4
      sample_runs = {}
      sample_fastqs.each do |file|
        machine, lane, sample = File.basename(file).split("_")
        sample_runs[[machine,lane]*"."] ||= []
        sample_runs[[machine,lane]*"."] << file
      end
      jobs = []
      num = 1
      sample_runs.each do |run_code,files|
        run_options = {}
        run_options[:read_group_name] = run_code + "." + sample_name
        run_options[:platform_unit] = run_code
        run_options[:fastq1] = files.sort.first
        run_options[:fastq2] = files.sort.last
        jobs << {:task => :BAM, :inputs => options.merge(run_options), :jobname => jobname + "_multiplex_" + num.to_s}
        num += 1
      end
      jobs
    else
      raise "Number of fastq is not 2 or 4: #{Misc.fingerprint sample_fastqs}"
    end
  end
  dep HTS, :BAM_multiplex, :compute => :ignore, :reference => Bellmunt.reference, :bam_files => :placeholder do |jobname,options,dependencies|
    if dependencies.flatten.length > 1
      {:jobname => jobname, :inputs => options.merge(:bam_files => dependencies.flatten.collect{|dep| dep.path})}
    else
      []
    end
  end
  dep_task :BAM, HTS, :BAM_rescore do |jobname,options, dependencies|
    if (mutiplex = dependencies.flatten.select{|dep| dep.task_name == :BAM_multiplex}.first)
      {:inputs => options.merge("HTS#BAM_duplicates" =>  mutiplex), :jobname => jobname + '_multiplexed'}
    else
      []
    end
  end

  dep :BAM
  extension :vcf
  dep_task :mutect2_snv, HTS, :mutect2_clean, :tumor => :BAM, :normal => nil, :reference => :placeholder,
    :interval_list => Bellmunt.interval_list,
    :pon => Bellmunt.pon,
		:germline_resource => Bellmunt.germline_resource,
    :af_not_in_resource => Bellmunt.af_not_in_resource do |jobname,options,dependencies|
      sample = jobname
      reference = dependencies.first.recursive_inputs[:reference]

      options[:reference] = Bellmunt.reference

      {:inputs => options, :jobname => jobname}
  end

  dep :BAM
  extension :vcf
  dep_task :strelka, HTS, :strelka, :tumor => :BAM, :normal => nil, :reference => :placeholder do |jobname,options,dependencies|
      sample = jobname
      reference = dependencies.first.recursive_inputs[:reference]

      options[:reference] = Bellmunt.reference

      {:inputs => options, :jobname => jobname}
  end

  dep :mutect2_snv, :compute => :produce
  dep_task :affected_genes, Sequence, :affected_genes, :mutations => :mutect2_snv, :organism => "Hsa/feb2014", :vcf => true

  dep :mutect2_snv, :compute => :produce
  dep_task :snv_annotations, Sample, :genomic_mutation_annotations, :file => :mutect2_snv, :vcf => true, :organism => Bellmunt.organism

  dep_task :OptiType, HLA, :OptiType do |jobname, options|
    sample_name = jobname

    options[:sample_name] = sample_name.gsub('_', '.')

    options[:reference] = Bellmunt.reference
    sample_fastqs = []
    Bellmunt.path.glob("*.fastq.gz").each do |file|
      basename = File.basename(file)
      sample_fastqs << file if basename.split("_")[2] == sample_name
    end
    {:inputs => options.merge(:files => sample_fastqs), :jobname => jobname}
  end

  dep :BAM
  dep_task :SOAPHLA, HLA, :SOAPHLA, :BAM => :BAM, :reference => Bellmunt.reference

  dep :mutect2_snv
  dep_task :genomic_mutations, Sequence, :genomic_mutations, :vcf_file => :mutect2_snv

  task :organism => :text do
    Bellmunt.organism
  end

end

module Sample

  dep Bellmunt, :mutect2_snv
  dep_task :genomic_mutations, Sequence, :genomic_mutations, :vcf_file => :mutect2_snv

  dep Bellmunt, :mutect2_snv
  dep_task :expanded_vcf, Sequence, :expanded_vcf, :vcf_file => :mutect2_snv

  extension :bam
  dep_task :BAM, Bellmunt, :BAM

  task :organism => :text do
    Bellmunt.organism
  end

  task :reference => :text do
    Bellmunt.reference
  end

  task :watson => :boolean do
    true
  end
end

#require 'Bellmunt/tasks/basic.rb'

#require 'rbbt/knowledge_base/Bellmunt'
#require 'rbbt/entity/Bellmunt'

