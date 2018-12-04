require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/Bellmunt'

Workflow.require_workflow "HTS"
Workflow.require_workflow "Sequence"
#Workflow.require_workflow "HLA"

module Bellmunt
  extend Workflow
  def self.reference
		'b37'
  end

  def self.path
    Rbbt.share.projects.Bellmunt
  end

  def self.interval_list
		Rbbt.share.projects.Bellmunt.data['Padded.hs37d5.bed'].find
  end

  def self.pon
		Rbbt.share.projects.Bellmunt.data['pon.vcf'].find
  end

  def self.af_not_in_resource
    0.0000025
  end

  def self.germline_resource
    GATK.known_sites.b37["af-only-gnomad.vcf"].produce.find
	end

  dep HTS, :BAM_rescore, :fastq1 => :placeholder, :fastq2 => :placeholder, :reference => :placeholder,
    :interval_list => Bellmunt.interval_list do |jobname,options|

    sample = jobname
    options[:fastq1] = self.path[sample + "_1.fastq.gz"].find
    options[:fastq2] = self.path[sample + "_2.fastq.gz"].find

    options[:reference] = Bellmunt.reference

    {:inputs => options, :jobname => jobname}
  end
  task :BAM => :binary do 
    Open.rm self.path
    Open.ln_s step(:BAM_rescore).path, self.path
    nil
  end

  dep :BAM
  dep HTS, :mutect2_clean, :tumor => :BAM, :normal => nil, :reference => :placeholder,
    :interval_list => Bellmunt.interval_list,
    :pon => Bellmunt.pon,
		:germline_resource => Bellmunt.germline_resource,
    :af_not_in_resource => Bellmunt.af_not_in_resource do |jobname,options,dependencies|
      sample = jobname
      reference = dependencies.first.recursive_inputs[:reference]

      options[:reference] = reference

      {:inputs => options, :jobname => jobname}
  end
  task :somatic_mutations => :tsv do
    TSV.get_stream step(:mutect2_clean)
  end

  dep :somatic_mutations
  dep Sequence, :affected_genes, :mutations => :somatic_mutations, :organism => "Hsa/feb2014"
  task :affected_genes => :array do
    TSV.get_stream step(:affected_genes)
  end

end

#require 'Bellmunt/tasks/basic.rb'

#require 'rbbt/knowledge_base/Bellmunt'
#require 'rbbt/entity/Bellmunt'

