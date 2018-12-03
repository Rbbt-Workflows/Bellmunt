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
    BWA.references.hs37d5['hs37d5.fa'].find
  end
  def self.path
    Rbbt.share.projects.Bellmunt
  end
  def self.interval_list
    Rbbt.share.projects.Bellmunt['intervals.bed'].find
  end

  dep HTS, :BAM_rescore do |jobname,options|

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
  dep HTS, :mutect2_clean, :tumor => :BAM do |jobname,options|
   sample = jobname
   options[:reference] = Bellmunt.reference
   options[:fastq1] = self.path[sample + "_1.fastq.gz"].find
   options[:fastq2] = self.path[sample + "_2.fastq.gz"].find
   {:inputs => options}
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

