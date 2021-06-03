cwlVersion: v1.0
class: Workflow
label: GATK Best Practice Data Pre-processing 4.1.0.0
doc: |-
  **BROAD Best Practice Data Pre-processing Workflow 4.1.0.0**  is used to prepare data for variant calling analysis. 

  It can be divided into two major segments: alignment to reference genome and data cleanup operations that correct technical biases [1].

  *A list of all inputs and parameters with corresponding descriptions can be found at the bottom of this page.*

  ***Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.***


  ### Common Use Cases

  * **BROAD Best Practice Data Pre-processing Workflow 4.1.0.0**  is designed to operate on individual samples.
  * Resulting BAM files are ready for variant calling analysis and can be further processed by other BROAD best practice pipelines, like **Generic germline short variant per-sample calling workflow** [2], **Somatic CNVs workflow** [3] and **Somatic SNVs+Indel workflow** [4].


  ### Changes Introduced by Seven Bridges

  This pipeline represents the CWL implementation of BROADs [original WDL file](https://github.com/gatk-workflows/gatk4-data-processing/pull/14) available on github. Minor differences are introduced in order to successfully adapt to the Seven Bridges Platform. These differences are listed below:
  * **SamToFastqAndBwaMem** step is divided into elementary steps: **SamToFastq** - converting unaligned BAM file to interleaved  FASTQ file, **BWA Mem** - performing alignment and **Samtools View** - used for converting SAM file to BAM.
  *  A boolean parameter **Ignore default RG ID** is added to **BWA MEM Bundle** tool. When used, this parameter ensures that **BWA MEM Bundle** does not add read group information (RG) in the BAM file. Instead, RG ID information obtained from uBAM is added by **GATK MergeBamAlignment** afterwards.* **SortAndFixTags** is divided into elementary steps: **SortSam** and **SetNmMdAndUqTags**
  * Added **SBG Lines to Interval List**: this tool is used to adapt results obtained with **CreateSequenceGroupingTSV**  for platform execution, more precisely for scattering.



  ### Common Issues and Important Notes

  * **BROAD Best Practice Data Pre-processing Workflow 4.1.0.0**  expects unmapped BAM file format as the main input.
  * **Input Alignments** (`--in_alignments`) - provided an unmapped BAM (uBAM) file should be in query-sorter order and all reads must have RG tags. Also, input uBAM files must pass validation by **ValidateSamFile**.
  * For each tool in the workflow, equivalent parameter settings to the one listed in the corresponding WDL file are set as defaults. 

  ### Performance Benchmarking
  Since this CWL implementation is meant to be equivalent to GATKs original WDL, there are no additional optimization steps beside instance and storage definition. 
  The c5.9xlarge AWS instance hint is used for WGS inputs and attached storage is set to 1.5TB.
  In the table given below one can find results of test runs for WGS and WES samples. All calculations are performed with reference files corresponding to assembly 38.

  *Cost can be significantly reduced by spot instance usage. Visit the [knowledge center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*

  | Input Size | Experimental Strategy | Coverage| Duration | Cost (spot) | AWS Instance Type |
  | --- | --- | --- | --- | --- | --- | 
  | 6.6 GiB | WES | 70 |1h 19min | $2.61 | c5.9 |
  |3.4 GiB | WES |  40 | 42min   | $1.40 | c5.9 |
  | 111.3 GiB| WGS | 30 |22h 41min | $43.86 | c5.9 |
  | 37.2 GiB  | WGS | 10 | 4h 21min | $14.21 | c5.9 |



  ### API Python Implementation
  The app's draft task can also be submitted via the **API**. In order to learn how to get your **Authentication token** and **API endpoint** for corresponding platform visit our [documentation](https://github.com/sbg/sevenbridges-python#authentication-and-configuration).

  ```python
  # Initialize the SBG Python API
  from sevenbridges import Api
  api = Api(token="enter_your_token", url="enter_api_endpoint")
  # Get project_id/app_id from your address bar. Example: https://igor.sbgenomics.com/u/your_username/project/app
  project_id = "your_username/project"
  app_id = "your_username/project/app"
  # Replace inputs with appropriate values
  inputs = {
  	"in_alignments": list(api.files.query(project=project_id, names=["<unaligned_bam>"])), 
  	"reference_index_tar": api.files.query(project=project_id, names=["Homo_sapiens_assembly38.fasta.tar"])[0], 
  	"in_reference": api.files.query(project=project_id, names=["Homo_sapiens_assembly38.fasta"])[0], 
  	"ref_dict": api.files.query(project=project_id, names=["Homo_sapiens_assembly38.dict"])[0],
  	"known_snps": api.files.query(project=project_id, names=["Homo_sapiens_assembly38.dbsnp.vcf"])[0],
          "known_sites": list(api.files.query(project=project_id, names=["Homo_sapiens_assembly38.known_indels.vcf", “Mills_and_1000G_gold_standard.indels.hg38.vcf”, “Homo_sapiens_assembly38.dbsnp.vcf”
  ]))}
  # Creates draft task
  task = api.tasks.create(name="BROAD Best Practice Data Pre-processing Workflow 4.1.0.0 - API Run", project=project_id, app=app_id, inputs=inputs, run=False)
  ```

  Instructions for installing and configuring the API Python client, are provided on [github](https://github.com/sbg/sevenbridges-python#installation). For more information about using the API Python client, consult [the client documentation](http://sevenbridges-python.readthedocs.io/en/latest/). **More examples** are available [here](https://github.com/sbg/okAPI).

  Additionally, [API R](https://github.com/sbg/sevenbridges-r) and [API Java](https://github.com/sbg/sevenbridges-java) clients are available. To learn more about using these API clients please refer to the [API R client documentation](https://sbg.github.io/sevenbridges-r/), and [API Java client documentation](https://docs.sevenbridges.com/docs/java-library-quickstart).


  ### References

  [1] [Data Pre-processing](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165)
  [2] [Generic germline short variant per-sample calling](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145)
  [3] [Somatic CNVs](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11147)
  [4] [Somatic SNVs+Indel pipeline ](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146)
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ScatterFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: in_alignments
  label: Input alignments
  doc: Input alignments files in unmapped BAM format.
  type: File[]
  sbg:fileTypes: SAM, BAM
  sbg:x: -648.1359252929688
  sbg:y: 25.01337432861328
- id: reference_index_tar
  label: BWA index archive
  doc: FASTA reference or BWA index archive.
  type: File
  sbg:fileTypes: TAR
  sbg:suggestedValue:
    name: GRCh38_primary_assembly_plus_ebv_alt_decoy_hla.fasta.tar
    class: File
    path: 5b6ace6e7550b4c330563856
  sbg:x: -583.3368530273438
  sbg:y: 259.1632995605469
- id: in_reference
  label: FASTA reference
  doc: Input reference in FASTA format.
  type: File
  secondaryFiles:
  - .fai
  - ^.dict
  sbg:fileTypes: FASTA, FA
  sbg:suggestedValue:
    name: Homo_sapiens_assembly38.fasta
    class: File
    path: 5772b6c7507c1752674486d1
  sbg:x: -447.3492126464844
  sbg:y: 555
- id: ref_dict
  label: DICT file
  doc: DICT file corresponding to the FASTA reference.
  type: File
  sbg:fileTypes: DICT
  sbg:suggestedValue:
    name: Homo_sapiens_assembly38.dict
    class: File
    path: 5c9ce4687369c402ac8a3c41
  sbg:x: 599.5844116210938
  sbg:y: -34.96286392211914
- id: known_sites
  label: Known sites
  doc: |-
    One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.  This argument must be specified at least once.
  type: File[]
  secondaryFiles:
  - |-
    ${
        var in_sites = self;
        if (in_sites.nameext == ".gz" || in_sites.nameext == '.GZ') {
                var tmp = in_sites.basename.slice(-7);
                if(tmp.toLowerCase() == '.vcf.gz') {
                    return in_sites.basename + ".tbi";  
                }
        }
        else if (in_sites.nameext == '.vcf' || in_sites.nameext == '.VCF' || in_sites.nameext == '.bed' || in_sites.nameext == '.BED') {
            return in_sites.basename + ".idx";
        }
        return in_sites.basename + ".idx";
    }
  sbg:fileTypes: VCF, VCF.GZ, BED
  sbg:x: 867.6756591796875
  sbg:y: 580.4737548828125

outputs:
- id: out_alignments
  label: Output BAM file
  doc: Output BAM file.
  type: File?
  outputSource:
  - gatk_gatherbamfiles_4_1_0_0/out_alignments
  sbg:fileTypes: BAM
  sbg:x: 2052.86767578125
  sbg:y: 289.4576416015625
- id: out_md5
  label: MD5 file
  doc: MD5 sum of the output BAM file.
  type: File?
  outputSource:
  - gatk_gatherbamfiles_4_1_0_0/out_md5
  sbg:fileTypes: MD5
  sbg:x: 2048
  sbg:y: 114.24113464355469
- id: out_duplication_metrics
  label: Duplication metrics
  doc: Duplication metrics file produced by GATK MarkDuplicates.
  type: File
  outputSource:
  - gatk_markduplicates_4_1_0_0/output_metrics
  sbg:fileTypes: METRICS
  sbg:x: 457.1893615722656
  sbg:y: -51.47343826293945

steps:
- id: gatk_markduplicates_4_1_0_0
  label: GATK MarkDuplicates
  in:
  - id: assume_sort_order
    default: queryname
  - id: in_alignments
    source:
    - gatk_mergebamalignment_4_1_0_0/out_alignments
  - id: optical_duplicate_pixel_distance
    default: 2500
  - id: validation_stringency
    default: SILENT
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK MarkDuplicates
    doc: |-
      The **GATK  MarkDuplicates** tool identifies duplicate reads in a BAM or SAM file.

      This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates [1].

      The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in the SAM/BAM file. The **Barcode tag** (`--BARCODE_TAG`) option is available to facilitate duplicate marking using molecular barcodes. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).


      ###Common Use Cases

      * The **GATK MarkDuplicates** tool requires the BAM or SAM file on its **Input BAM/SAM file** (`--INPUT`) input. The tool generates a new SAM or BAM file on its **Output BAM/SAM** output, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, please see the following [blog post](https://software.broadinstitute.org/gatk/blog?id=7019) for additional information. **MarkDuplicates** also produces a metrics file on its **Output metrics file** output, indicating the numbers of duplicates for both single and paired end reads.

      * The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.

      * If desired, duplicates can be removed using the **Remove duplicates** (`--REMOVE_DUPLICATES`) and **Remove sequencing duplicates** ( `--REMOVE_SEQUENCING_DUPLICATES`) options.

      * Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type of duplicate. To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output of a SAM/BAM file. Invoking the **Tagging policy** ( `--TAGGING_POLICY`) option, you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the output SAM/BAM file will have values for the 'DT' tag (depending on the invoked **TAGGING_POLICY** option), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). 

      * This tool uses the **Read name regex** (`--READ_NAME_REGEX`) and the **Optical duplicate pixel distance** (`--OPTICAL_DUPLICATE_PIXEL_DISTANCE`) options as the primary methods to identify and differentiate duplicate types. Set **READ_NAME_REGEX** to null to skip optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are extremely large and estimating library complexity is not an aim. Note that without optical duplicate counts, library size estimation will be inaccurate.

      * Usage example:

      ```
      gatk MarkDuplicates \
            --INPUT input.bam \
            --OUTPUT marked_duplicates.bam \
            --METRICS_FILE marked_dup_metrics.txt
      ```

      ###Changes Introduced by Seven Bridges

      * All output files will be prefixed using the **Output prefix** parameter. In case **Output prefix** is not provided, output prefix will be the same as the Sample ID metadata from the **Input SAM/BAM file**, if the Sample ID metadata exists. Otherwise, output prefix will be inferred from the **Input SAM/BAM** filename. This way, having identical names of the output files between runs is avoided. Moreover,  **dedupped** will be added before the extension of the output file name. 

      * The user has a possibility to specify the output file format using the **Output file format** option. Otherwise, the output file format will be the same as the format of the input file.

      ###Common Issues and Important Notes

      * None

      ###Performance Benchmarking

      Below is a table describing runtimes and task costs of **GATK MarkDuplicates** for a couple of different samples, executed on the AWS cloud instances:

      | Experiment type |  Input size | Duration |  Cost | Instance (AWS) | 
      |:--------------:|:------------:|:--------:|:-------:|:---------:|
      |     RNA-Seq     |  1.8 GB |   3min   | ~0.02$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     |  5.3 GB |   9min   | ~0.06$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 8.8 GB |  16min  | ~0.11$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 17 GB |  30min  | ~0.20$ | c4.2xlarge (8 CPUs) |

      *Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*

      ###References

      [1] [GATK MarkDuplicates](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_markduplicates_MarkDuplicates.php)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: "${\n    return inputs.cpu_per_job ? inputs.cpu_per_job : 1;\n}"
      ramMin: |-
        ${
            var memory = 4096;
            if (inputs.memory_per_job) 
            {
                memory = inputs.memory_per_job;
            }
            if (inputs.memory_overhead_per_job)
            {
                memory += inputs.memory_overhead_per_job;
            }
            return memory;
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

    inputs:
    - id: add_pg_tag_to_reads
      label: Add PG tag to reads
      doc: Add PG tag to each read in a SAM or BAM file.
      type:
      - 'null'
      - name: add_pg_tag_to_reads
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --ADD_PG_TAG_TO_READS
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: assume_sort_order
      label: Assume sort order
      doc: |-
        If not null, assume that the input file has this order even if the header says otherwise. Cannot be used in conjuction with argument(s) ASSUME_SORTED (AS).
      type:
      - 'null'
      - name: assume_sort_order
        type: enum
        symbols:
        - unsorted
        - queryname
        - coordinate
        - duplicate
        - unknown
      inputBinding:
        prefix: --ASSUME_SORT_ORDER
        position: 4
        shellQuote: false
      sbg:altPrefix: -ASO
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: assume_sorted
      label: Assume sorted
      doc: |-
        If true, assume that the input file is coordinate sorted even if the header says otherwise. Deprecated, used ASSUME_SORT_ORDER=coordinate instead. Exclusion: This argument cannot be used at the same time as ASSUME_SORT_ORDER (ASO).
      type: boolean?
      inputBinding:
        prefix: --ASSUME_SORTED
        position: 4
        shellQuote: false
      sbg:altPrefix: -AS
      sbg:category: Optional arguments
      sbg:toolDefaultValue: 'false'
    - id: barcode_tag
      label: Barcode tag
      doc: Barcode SAM tag (ex. BC for 10x genomics).
      type: string?
      inputBinding:
        prefix: --BARCODE_TAG
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: clear_dt
      label: Clear DT
      doc: |-
        Clear DT tag from input SAM records. Should be set to false if input SAM doesn't have this tag.
      type:
      - 'null'
      - name: clear_dt
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --CLEAR_DT
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: comment
      label: Comment
      doc: Comment(s) to include in the output file's header.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--COMMENT', self[i]);
                      
                  }
                  return cmd.join(' ');
              }
          }
        shellQuote: false
      sbg:altPrefix: -CO
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: compression_level
      label: Compression level
      doc: Compression level for all compressed files created (e.g. BAM and VCF).
      type: int?
      inputBinding:
        prefix: --COMPRESSION_LEVEL
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
    - id: create_index
      label: Create index
      doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
      type: boolean?
      inputBinding:
        prefix: --CREATE_INDEX
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: duplex_umi
      label: Duplex UMI
      doc: |-
        Treat UMIs as being duplex stranded. This option requires that the UMI consist of two equal length strings that are separated by a hyphen (e.g. 'ATC-GTC'). Reads are considered duplicates if, in addition to standard definition, have identical normalized UMIs. A UMI from the 'bottom' strand is normalized by swapping its content around the hyphen (eg. ATC-GTC becomes GTC-ATC). A UMI from the 'top' strand is already normalized as it is. Both reads from a read pair considered top strand if the read 1 unclipped 5' coordinate is less than the read 2 unclipped 5' coordinate. All chimeric reads and read fragments are treated as having come from the top strand. With this option it is required that the BARCODE_TAG hold non-normalized UMIs.
      type: boolean?
      inputBinding:
        prefix: --DUPLEX_UMI
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: duplicate_scoring_strategy
      label: Duplicate scoring strategy
      doc: The scoring strategy for choosing the non-duplicate among candidates.
      type:
      - 'null'
      - name: duplicate_scoring_strategy
        type: enum
        symbols:
        - SUM_OF_BASE_QUALITIES
        - TOTAL_MAPPED_REFERENCE_LENGTH
        - RANDOM
      inputBinding:
        prefix: --DUPLICATE_SCORING_STRATEGY
        position: 4
        shellQuote: false
      sbg:altPrefix: -DS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: SUM_OF_BASE_QUALITIES
    - id: in_alignments
      label: Input BAM/SAM file
      doc: Input SAM or BAM files to analyze. Must be coordinate sorted.
      type: File[]
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              var in_files = [].concat(inputs.in_alignments);
              if (in_files)
              {
                  var cmd = [];
                  for (var i = 0; i < in_files.length; i++) 
                  {
                      cmd.push('--INPUT', in_files[i].path);
                  }
                  return cmd.join(' ');
              }
          }
        shellQuote: false
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM
    - id: max_file_handles_for_read_ends_map
      label: Max file handles for read ends map
      doc: |-
        Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a unix system.
      type: int?
      inputBinding:
        prefix: --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP
        position: 4
        shellQuote: false
      sbg:altPrefix: -MAX_FILE_HANDLES
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '8000'
    - id: max_optical_duplicate_set_size
      label: Max optical duplicate set size
      doc: |-
        This number is the maximum size of a set of duplicate reads for which we will attempt to determine which are optical duplicates. Please be aware that if you raise this value too high and do encounter a very large set of duplicate reads, it will severely affect the runtime of this tool. To completely disable this check, set the value to -1.
      type: int?
      inputBinding:
        prefix: --MAX_OPTICAL_DUPLICATE_SET_SIZE
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '300000'
    - id: max_records_in_ram
      label: Max records in RAM
      doc: |-
        When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
      type: int?
      inputBinding:
        prefix: --MAX_RECORDS_IN_RAM
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
    - id: memory_overhead_per_job
      label: Memory overhead per job
      doc: |-
        This input allows a user to set the desired overhead memory when running a tool or adding it to a workflow. This amount will be added to the Memory per job in the Memory requirements section but it will not be added to the -Xmx parameter leaving some memory not occupied which can be used as stack memory (-Xmx parameter defines heap memory). This input should be defined in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      type: int?
      sbg:category: Platform Options
    - id: memory_per_job
      label: Memory per job
      doc: |-
        This input allows a user to set the desired memory requirement when running a tool or adding it to a workflow. This value should be propagated to the -Xmx parameter too.This input should be defined in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      type: int?
      sbg:category: Platform Options
    - id: molecular_identifier_tag
      label: Molecular identifier tag
      doc: |-
        SAM tag to uniquely identify the molecule from which a read was derived. Use of this option requires that the BARCODE_TAG option be set to a non null value.
      type: string?
      inputBinding:
        prefix: --MOLECULAR_IDENTIFIER_TAG
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: optical_duplicate_pixel_distance
      label: Optical duplicate pixel distance
      doc: |-
        The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the illumina platform. For the patterned flowcell models, 2500 is moreappropriate. For other platforms and models, users should experiment to find what works best.
      type: int?
      inputBinding:
        prefix: --OPTICAL_DUPLICATE_PIXEL_DISTANCE
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '100'
    - id: program_group_command_line
      label: Program group command line
      doc: |-
        Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically.
      type: string?
      inputBinding:
        prefix: --PROGRAM_GROUP_COMMAND_LINE
        position: 4
        shellQuote: false
      sbg:altPrefix: -PG_COMMAND
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: program_group_name
      label: Program group name
      doc: Value of PN tag of PG record to be created.
      type: string?
      inputBinding:
        prefix: --PROGRAM_GROUP_NAME
        position: 4
        shellQuote: false
      sbg:altPrefix: -PG_NAME
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: MarkDuplicates
    - id: program_group_version
      label: Program group version
      doc: |-
        Value of VN tag of PG record to be created. If not specified, the version will be detected automatically.
      type: string?
      inputBinding:
        prefix: --PROGRAM_GROUP_VERSION
        position: 4
        shellQuote: false
      sbg:altPrefix: -PG_VERSION
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: program_record_id
      label: Program record id
      doc: |-
        The program record ID for the @PG record(s) created by this program. Set to null to disable PG record creation.  This string may have a suffix appended to avoid collision with other program record IDs.
      type: string?
      inputBinding:
        prefix: --PROGRAM_RECORD_ID
        position: 4
        shellQuote: false
      sbg:altPrefix: -PG
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: MarkDuplicates
    - id: read_name_regex
      label: Read name regex
      doc: |-
        MarkDuplicates can use the tile and cluster positions to estimate the rate of optical duplication in addition to the dominant source of duplication, PCR, to provide a more accurate estimation of library size. By default (with no READ_NAME_REGEX specified), MarkDuplicates will attempt to extract coordinates using a split on ':' (see note below). Set READ_NAME_REGEX to 'null' to disable optical duplicate detection. Note that without optical duplicate counts, library size estimation will be less accurate. If the read name does not follow a standard illumina colon-separation convention, but does contain tile and x,y coordinates, a regular expression can be specified to extract three variables: tile/region, x coordinate and y coordinate from a read name. The regular expression must contain three capture groups for the three variables, in order. It must match the entire read name. e.g. if field names were separated by semi-colon (';') this example regex could be specified (?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$ Note that if no READ_NAME_REGEX is specified, the read name is split on ':'. For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values.
      type: string?
      inputBinding:
        prefix: --READ_NAME_REGEX
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
    - id: read_one_barcode_tag
      label: Read one barcode tag
      doc: Read one barcode SAM tag (ex. BX for 10x Genomics).
      type: string?
      inputBinding:
        prefix: --READ_ONE_BARCODE_TAG
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: read_two_barcode_tag
      label: Read two barcode tag
      doc: Read two barcode SAM tag (ex. BX for 10x Genomics).
      type: string?
      inputBinding:
        prefix: --READ_TWO_BARCODE_TAG
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: remove_duplicates
      label: Remove duplicates
      doc: |-
        If true do not write duplicates to the output file instead of writing them with appropriate flags set.
      type: boolean?
      inputBinding:
        prefix: --REMOVE_DUPLICATES
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: remove_sequencing_duplicates
      label: Remove sequencing duplicates
      doc: |-
        If true remove 'optical' duplicates and other duplicates that appear to have arisen from the sequencing process instead of the library preparation process, even if REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored.
      type: boolean?
      inputBinding:
        prefix: --REMOVE_SEQUENCING_DUPLICATES
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: sorting_collection_size_ratio
      label: Sorting collection size ratio
      doc: |-
        This number, plus the maximum RAM available to the JVM, determine the memory footprint used by some of the sorting collections. If you are running out of memory, try reducing this number.
      type: float?
      inputBinding:
        prefix: --SORTING_COLLECTION_SIZE_RATIO
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0.25'
    - id: tag_duplicate_set_members
      label: Tag duplicate set members
      doc: |-
        If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG (DS), indicates the size of the duplicate set. The smallest possible DS value is 2 which occurs when two reads map to the same portion of the reference only one of which is marked as duplicate. The second tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier for the duplicate set to which the record belongs. This identifier is the index-in-file of the representative read that was selected out of the duplicate set.
      type: boolean?
      inputBinding:
        prefix: --TAG_DUPLICATE_SET_MEMBERS
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: tagging_policy
      label: Tagging policy
      doc: Determines how duplicate types are recorded in the DT optional attribute.
      type:
      - 'null'
      - name: tagging_policy
        type: enum
        symbols:
        - DontTag
        - OpticalOnly
        - All
      inputBinding:
        prefix: --TAGGING_POLICY
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: DontTag
    - id: validation_stringency
      label: Validation stringency
      doc: |-
        Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
      type:
      - 'null'
      - name: validation_stringency
        type: enum
        symbols:
        - STRICT
        - LENIENT
        - SILENT
      inputBinding:
        prefix: --VALIDATION_STRINGENCY
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
    - id: output_prefix
      label: Output prefix
      doc: Output file name prefix.
      type: string?
      sbg:category: Optional Arguments
    - id: output_file_format
      label: Output file format
      doc: Output file format
      type:
      - 'null'
      - name: output_file_format
        type: enum
        symbols:
        - bam
        - sam
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: BAM
    - id: cpu_per_job
      label: CPU per job
      doc: |-
        This input allows a user to set the desired CPU requirement when running a tool or adding it to a workflow.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'

    outputs:
    - id: out_alignments
      label: Output BAM/SAM file
      doc: Output BAM/SAM file which contains marked records.
      type: File?
      secondaryFiles:
      - |-
        ${ 
           if (inputs.create_index)   {
               return [self.basename + ".bai", self.nameroot + ".bai"]
           }  else {
               return []; 
          }
        }
      outputBinding:
        glob: '*am'
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: BAM, SAM
    - id: output_metrics
      label: Output metrics file
      doc: Output duplication metrics file.
      type: File
      outputBinding:
        glob: '*metrics'
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: METRICS

    baseCommand: []
    arguments:
    - prefix: ''
      position: 0
      valueFrom: /opt/gatk
      shellQuote: false
    - prefix: ''
      position: 1
      valueFrom: |-
        ${
            if (inputs.memory_per_job)
            {
                return "--java-options";
            }
            else {
                return ''; 
            }
        }
            
      shellQuote: false
    - prefix: ''
      position: 2
      valueFrom: |-
        ${
            if (inputs.memory_per_job) {
                return '\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\"';
            }
            else {
                return ''; 
            }
        }
      shellQuote: false
    - position: 3
      valueFrom: MarkDuplicates
      shellQuote: false
    - prefix: ''
      position: 4
      valueFrom: |-
        ${
            var in_alignments = [].concat(inputs.in_alignments);
            var output_ext = inputs.output_file_format ? "." + inputs.output_file_format : in_alignments[0].nameext;
            var output_prefix = '';
            if (inputs.output_prefix)
            {
                output_prefix = inputs.output_prefix;
            }
            else
            {
                if (in_alignments[0].metadata && in_alignments[0].metadata.sample_id)
                {
                    output_prefix = in_alignments[0].metadata.sample_id;
                }
                else
                {
                    output_prefix = in_alignments[0].nameroot.split('.')[0];
                }
            }
            return "--OUTPUT " + output_prefix + ".dedupped" + output_ext;
        }
      shellQuote: false
    - prefix: ''
      position: 4
      valueFrom: |-
        ${
            var in_alignments = [].concat(inputs.in_alignments);
            var output_prefix = '';  

            if (inputs.output_prefix)
            {
                output_prefix = inputs.output_prefix;
            }
            else
            {
                if (in_alignments[0].metadata && in_alignments[0].metadata.sample_id)
                {
                    output_prefix = in_alignments[0].metadata.sample_id;
                }
                else
                {
                    output_prefix = in_alignments[0].nameroot.split('.')[0];
                }
            }
            return "--METRICS_FILE " + output_prefix + ".dedupped.metrics";
        }
      shellQuote: false
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-markduplicates-4-1-0-0/12
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    sbg:content_hash: a112438cd40b078b2fbf816496a7cabec5688e19c781aac7f79a1de917e0eabfb
    sbg:contributors:
    - uros_sipetic
    - nemanja.vucic
    - veliborka_josipovic
    - nens
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/26
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552668097
    sbg:id: h-b4115186/h-46b74571/h-420c567b/0
    sbg:image_url:
    sbg:latestRevision: 12
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - id: https://software.broadinstitute.org/gatk/
      label: Homepage
    - id: https://github.com/broadinstitute/gatk/
      label: Source Code
    - id: |-
        https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
      label: Download
    - id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
      label: Publications
    - id: |-
        https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_markduplicates_MarkDuplicates.php
      label: Documentation
    sbg:modifiedBy: uros_sipetic
    sbg:modifiedOn: 1562416183
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 12
    sbg:revisionNotes: |-
      Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/26
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1552668097
      sbg:revision: 0
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/9
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554492835
      sbg:revision: 1
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/13
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554720881
      sbg:revision: 2
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/14
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554999255
      sbg:revision: 3
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/15
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1555945044
      sbg:revision: 4
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/17
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557734534
      sbg:revision: 5
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/18
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000580
      sbg:revision: 6
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/19
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351536
      sbg:revision: 7
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/21
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558447931
      sbg:revision: 8
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/22
    - sbg:modifiedBy: nemanja.vucic
      sbg:modifiedOn: 1559750423
      sbg:revision: 9
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/23
    - sbg:modifiedBy: nemanja.vucic
      sbg:modifiedOn: 1559751034
      sbg:revision: 10
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/24
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1561632463
      sbg:revision: 11
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/25
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1562416183
      sbg:revision: 12
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/26
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: out_alignments
  - id: output_metrics
  sbg:x: 252.3874969482422
  sbg:y: 88.93749237060547
- id: bwa_mem_bundle_0_7_15
  label: BWA MEM Bundle
  in:
  - id: verbose_level
    default: '3'
  - id: smart_pairing_in_input_fastq
    default: true
  - id: input_reads
    source:
    - gatk_samtofastq_4_1_0_0/out_reads
  - id: num_input_bases_in_each_batch
    default: 100000000
  - id: use_soft_clipping
    default: true
  - id: threads
    default: 16
  - id: output_header
    default: false
  - id: reference_index_tar
    source: reference_index_tar
  - id: output_format
    default: SAM
  - id: mapQ_of_suplementary
    default: false
  - id: ignore_default_rg_id
    default: true
  scatter:
  - input_reads
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: BWA MEM Bundle 0.7.15 CWL1.0
    doc: |-
      BWA-MEM is an algorithm designed for aligning sequence reads onto a large reference genome. BWA-MEM is implemented as a component of BWA. The algorithm can automatically choose between performing end-to-end and local alignments. BWA-MEM is capable of outputting multiple alignments, and finding chimeric reads. It can be applied to a wide range of read lengths, from 70 bp to several megabases. 

      *A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*


      ## Common Use Cases
      In order to obtain possibilities for additional fast processing of aligned reads, **Biobambam2 sortmadup** (2.0.87) tool is embedded together into the same package with BWA-MEM (0.7.15).

      In order to obtain possibilities for additional fast processing of aligned reads, **Biobambam2** (2.0.87) is embedded together with the BWA 0.7.15 toolkit into the **BWA-MEM Bundle 0.7.15 CWL1.0**.  Two tools are used (**bamsort** and **bamsormadup**) to allow the selection of three output formats (SAM, BAM, or CRAM), different modes of sorting (Quarryname/Coordinate sorting), and Marking/Removing duplicates that can arise during sample preparation e.g. library construction using PCR. This is done by setting the **Output format** and **PCR duplicate detection** parameters.
      - Additional notes:
          - The default **Output format** is coordinate sorted BAM (option **BAM**).
          - SAM and BAM options are query name sorted, while CRAM format is not advisable for data sorted by query name.
          - Coordinate Sorted BAM file in all options and CRAM Coordinate sorted output with Marked Duplicates come with the accompanying index file. The generated index name will be the same as the output alignments file, with the extension BAM.BAI or CRAM.CRAI. However, when selecting the CRAM Coordinate sorted and CRAM Coordinate sorted output with Removed Duplicates, the generated files will not have the index file generated. This is a result of the usage of different Biobambam2 tools - **bamsort** does not have the ability to write CRAI files (only supports outputting BAI index files), while **bamsormadup** can write CRAI files.
          - Passing data from BWA-MEM to Biobambam2 tools has been done through the Linux piping which saves processing times (up to an hour of the execution time for whole-genome sample) of reading and writing of aligned reads into the hard drive. 
          - **BWA-MEM Bundle 0.7.15 CWL1** first needs to construct the FM-index  (Full-text index in Minute space) for the reference genome using the **BWA INDEX 0.7.17 CWL1.0** tool. The two BWA versions are compatible.

      ### Changes Introduced by Seven Bridges

      - **Aligned SAM/BAM/CRAM** file will be prefixed using the **Output SAM/BAM/CRAM file name** parameter. In case **Output SAM/BAM/CRAM file name** is not provided, the output prefix will be the same as the **Sample ID** metadata field from the file if the **Sample ID** metadata field exists. Otherwise, the output prefix will be inferred from the **Input reads** file names.
      -  The **Platform** metadata field for the output alignments will be automatically set to "Illumina" unless it is present in **Input reads** metadata, or given through **Read group header** or **Platform** input parameters. This will prevent possible errors in downstream analysis using the GATK toolkit.
      - If the **Read group ID** parameter is not defined, by default it will be set to ‘1’. If the tool is scattered within a workflow it will assign the **Read Group ID** according to the order of the scattered folders. This ensures a unique **Read Group ID** when processing multi-read group input data from one sample.

      ### Common Issues and Important Notes 
       
      - For input reads FASTQ files of total size less than 10 GB we suggest using the default setting for parameter **Total memory** of 15GB, for larger files we suggest using 58 GB of memory and 32 CPU cores.
      - When the desired output is a CRAM file without deduplication of the PCR duplicates, it is necessary to provide the FASTA Index file (FAI) as input.
      - Human reference genome version 38 comes with ALT contigs, a collection of diverged alleles present in some humans but not the others. Making effective use of these contigs will help to reduce mapping artifacts, however, to facilitate mapping these ALT contigs to the primary assembly, GRC decided to add to each contig long flanking sequences almost identical to the primary assembly. As a result, a naive mapping against GRCh38+ALT will lead to many mapQ-zero mappings in these flanking regions. Please use post-processing steps to fix these alignments or implement [steps](https://sourceforge.net/p/bio-bwa/mailman/message/32845712/) described by the author of the BWA toolkit.  
      - Inputs **Read group header** and **Insert string to header** need to be given in the correct format - under single-quotes.
      - BWA-MEM is not a splice aware aligner, so it is not the appropriate tool for mapping RNAseq to the genome. For RNAseq reads **Bowtie2 Aligner** and **STAR** are recommended tools. 
      - Input paired reads need to have the identical read names - if not, the tool will throw a ``[mem_sam_pe] paired reads have different names`` error.
      - This wrapper was tested and is fully compatible with cwltool v3.0.

      ### Performance Benchmarking

      Below is a table describing the runtimes and task costs on on-demand instances for a set of samples with different file sizes :

      | Input reads       | Size [GB] | Output format | Instance (AWS)           | Duration  | Cost   | Threads |
      |-------------------|-----------|---------------|--------------------------|-----------|--------|---------|
      | HG001-NA12878-30x | 2 x 23.8  | SAM           | c5.9xlarge (36CPU, 72GB) | 5h 12min  | $7.82  | 36      |
      | HG001-NA12878-30x | 2 x 23.8  | BAM           | c5.9xlarge (36CPU, 72GB) | 5h 16min  | $8.06  | 36      |
      | HG002-NA24385-50x | 2 x 66.4  | SAM           | c5.9xlarge (36CPU, 72GB) | 8h 33min  | $13.08 | 36      |


      *Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: |-
        ${
            var reads_size = 0
            // Calculate suggested number of CPUs depending of the input reads size
            if (inputs.input_reads.constructor == Array) {
                if (inputs.input_reads[1]) reads_size = inputs.input_reads[0].size + inputs.input_reads[1].size;
                else reads_size = inputs.input_reads[0].size;
            } else reads_size = inputs.input_reads.size;
            
            if (!reads_size) reads_size = 0;
            
            var GB_1 = 1024 * 1024 * 1024;
            var suggested_cpus = 0;
            if (reads_size < GB_1) suggested_cpus = 1;
            else if (reads_size < 10 * GB_1) suggested_cpus = 8;
            else suggested_cpus = 31;
            
            if (inputs.reserved_threads) return inputs.reserved_threads;
            else if (inputs.threads) return inputs.threads;
            else if (inputs.sambamba_threads) return inputs.sambamba_threads;
            else return suggested_cpus;
            
        }
      ramMin: |-
        ${
            var reads_size =0;
            // Calculate suggested number of CPUs depending of the input reads size
            if (inputs.input_reads.constructor == Array) {
                if (inputs.input_reads[1]) reads_size = inputs.input_reads[0].size + inputs.input_reads[1].size;
                else reads_size = inputs.input_reads[0].size;
            } else reads_size = inputs.input_reads.size;
            if (!reads_size) reads_size = 0;

            var GB_1 = 1024 * 1024 * 1024;
            var  suggested_memory = 0;
            if (reads_size < GB_1) suggested_memory = 4;
            else if (reads_size < 10 * GB_1) suggested_memory = 15;
            else suggested_memory = 58;
            
            if (inputs.total_memory) return inputs.total_memory * 1024;
            else if (inputs.sort_memory) return inputs.sort_memory * 1024;
            else return suggested_memory * 1024;
            
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/nens/bwa-0-7-15:0
    - class: InitialWorkDirRequirement
      listing:
      - $(inputs.reference_index_tar)
      - $(inputs.input_reads)
      - $(inputs.fasta_index)
    - class: InlineJavascriptRequirement
      expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


        var setMetadata = function(file, metadata) {
            if (!('metadata' in file)) {
                file['metadata'] = {}
            }
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };

    inputs:
    - id: drop_chains_fraction
      label: Drop chains fraction
      doc: |-
        Drop chains shorter than a given fraction (FLOAT) of the longest overlapping chain.
      type: float?
      inputBinding:
        prefix: -D
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '0.50'
    - id: verbose_level
      label: Verbose level
      doc: 'Select verbose level: 1=error, 2=warning, 3=message, 4+=debugging.'
      type:
      - 'null'
      - name: verbose_level
        type: enum
        symbols:
        - '1'
        - '2'
        - '3'
        - '4'
      inputBinding:
        prefix: -v
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
      sbg:toolDefaultValue: '3'
    - id: sort_memory
      label: Memory for BAM sorting
      doc: |-
        Amount of RAM [Gb] to give to the sorting algorithm (if not provided will be set to one-third of the total memory).
      type: int?
      sbg:category: Execution
    - id: wgs_hg38_mode_threads
      label: Optimize threads for HG38
      doc: Lower the number of threads if HG38 reference genome is used.
      type: int?
      sbg:category: Execution
      sbg:toolDefaultValue: 'False'
    - id: band_width
      label: Band width
      doc: Band width for banded alignment.
      type: int?
      inputBinding:
        prefix: -w
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '100'
    - id: smart_pairing_in_input_fastq
      label: Smart pairing in input FASTQ file
      doc: Smart pairing in input FASTQ file (ignoring in2.fq).
      type: boolean?
      inputBinding:
        prefix: -p
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
    - id: rg_library_id
      label: Library ID
      doc: |-
        Specify the identifier for the sequencing library preparation, which will be placed in RG line.
      type: string?
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Inferred from metadata
    - id: mate_rescue_rounds
      label: Mate rescue rounds
      doc: |-
        Perform at the most a given number (INT) of rounds of mate rescues for each read.
      type: string?
      inputBinding:
        prefix: -m
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '50'
    - id: reserved_threads
      label: Reserved number of threads on the instance
      doc: Reserved number of threads on the instance used by scheduler.
      type: int?
      sbg:category: Configuration
      sbg:toolDefaultValue: '1'
    - id: input_reads
      label: Input reads
      doc: Input sequence reads.
      type: File[]
      inputBinding:
        position: 105
        valueFrom: |-
          ${
              /// Set input reads in the correct order depending of the paired end from metadata

              // Set output file name
              function flatten(files){
                  var a = [];
                  for(var i=0;i<files.length;i++){
                      if(files[i]){
                          if(files[i].constructor == Array) a = a.concat(flatten(files[i]));
                          else a = a.concat(files[i]);}}
                  var b = a.filter(function (el) {return el != null;})
                  return b;}
              var files1 = [].concat(inputs.input_reads);
              var in_reads=flatten(files1);

              // Read metadata for input reads
              var read_metadata = in_reads[0].metadata;
              if (!read_metadata) read_metadata = [];

              var order = 0; // Consider this as normal order given at input: pe1 pe2

              // Check if paired end 1 corresponds to the first given read
              if (read_metadata == []) order = 0;
              else if ('paired_end' in read_metadata) {
                  var pe1 = read_metadata.paired_end;
                  if (pe1 != 1) order = 1; // change order
              }

              // Return reads in the correct order
              if (in_reads.length == 1) return in_reads[0].path; // Only one read present
              else if (in_reads.length == 2) {
                  if (order == 0) return in_reads[0].path + ' ' + in_reads[1].path;
                  else return in_reads[1].path + ' ' + in_reads[0].path;
              }
          }
        shellQuote: false
      sbg:category: Input files
      sbg:fileTypes: FASTQ, FASTQ.GZ, FQ, FQ.GZ
    - id: unpaired_read_penalty
      label: Unpaired read penalty
      doc: Penalty for an unpaired read pair.
      type: int?
      inputBinding:
        prefix: -U
        position: 4
        shellQuote: false
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '17'
    - id: clipping_penalty
      label: Clipping penalty
      doc: Penalty for 5'- and 3'-end clipping.
      type: int[]?
      inputBinding:
        prefix: -L
        position: 4
        separate: false
        itemSeparator: ','
        shellQuote: false
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '[5,5]'
    - id: select_seeds
      label: Select seeds
      doc: Look for internal seeds inside a seed longer than {-k} * FLOAT.
      type: float?
      inputBinding:
        prefix: -r
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '1.5'
    - id: score_for_a_sequence_match
      label: Score for a sequence match
      doc: Score for a sequence match, which scales options -TdBOELU unless overridden.
      type: int?
      inputBinding:
        prefix: -A
        position: 4
        shellQuote: false
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '1'
    - id: dropoff
      label: Dropoff
      doc: Off-diagonal X-dropoff.
      type: int?
      inputBinding:
        prefix: -d
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '100'
    - id: num_input_bases_in_each_batch
      label: Number of input bases to process
      doc: |-
        Process a given number (INT) of input bases in each batch regardless of nThreads (for reproducibility).
      type: int?
      inputBinding:
        prefix: -K
        position: 4
        shellQuote: false
    - id: total_memory
      label: Total memory
      doc: |-
        Total memory to be used by the tool in GB. It's the sum of BWA and BIOBAMBAM2 processes. For FASTQ files of a total size less than 10GB, we suggest using the default setting of 15GB, for larger files, we suggest using 58GB of memory (and 32CPU cores).
      type: int?
      sbg:category: Execution
      sbg:toolDefaultValue: '15'
    - id: gap_extension_penalties
      label: Gap extension
      doc: |-
        Gap extension penalty; a gap of size k cost '{-O} + {-E}*k'. 
        This array can't have more than two values.
      type: int[]?
      inputBinding:
        prefix: -E
        position: 4
        separate: false
        itemSeparator: ','
        shellQuote: false
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '[1,1]'
    - id: deduplication
      label: PCR duplicate detection
      doc: Use Biobambam2 for finding duplicates on sequence reads.
      type:
      - 'null'
      - name: deduplication
        type: enum
        symbols:
        - None
        - MarkDuplicates
        - RemoveDuplicates
      sbg:category: Biobambam2 parameters
      sbg:toolDefaultValue: MarkDuplicates
    - id: ignore_alt_file
      label: Ignore ALT file
      doc: |-
        Treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file).
      type: boolean?
      inputBinding:
        prefix: -j
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
    - id: rg_id
      label: Read group ID
      doc: Set read group ID.
      type: string?
      sbg:category: Configuration
      sbg:toolDefaultValue: '1'
    - id: use_soft_clipping
      label: Use soft clipping
      doc: Use soft clipping for supplementary alignments.
      type: boolean?
      inputBinding:
        prefix: -Y
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
    - id: output_in_xa
      label: Output in XA
      doc: |-
        If there are < number (INT) of hits with a score >80% of the max score, output all in XA. 
        This array should have no more than two values.
      type: int[]?
      inputBinding:
        prefix: -h
        position: 4
        separate: false
        itemSeparator: ','
        shellQuote: false
      sbg:category: BWA Input/output options
      sbg:toolDefaultValue: '[5, 200]'
    - id: rg_platform
      label: Platform
      doc: |-
        Specify the version of the technology that was used for sequencing, which will be placed in RG line.
      type:
      - 'null'
      - name: rg_platform
        type: enum
        symbols:
        - '454'
        - Helicos
        - Illumina
        - Solid
        - IonTorrent
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Inferred from metadata
    - id: threads
      label: Threads
      doc: |-
        The number of threads for BWA and Biobambam2 sort processes (both will use the given number).
      type: int?
      sbg:category: Execution
      sbg:toolDefaultValue: '8'
    - id: skip_pairing
      label: Skip pairing
      doc: Skip pairing; mate rescue performed unless -S also in use.
      type: boolean?
      inputBinding:
        prefix: -P
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
    - id: insert_string_to_header
      label: Insert string to header
      doc: Insert STR to output header if it starts with "@".
      type: string?
      inputBinding:
        prefix: -H
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
    - id: output_header
      label: Output header
      doc: Output the reference FASTA header in the XR tag.
      type: boolean?
      inputBinding:
        prefix: -V
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
    - id: seed_occurrence_for_the_3rd_round
      label: Seed occurrence
      doc: Seed occurrence for the 3rd round seeding.
      type: int?
      inputBinding:
        prefix: -y
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '20'
    - id: read_type
      label: Sequencing technology-specific settings
      doc: |-
        Sequencing technology-specific settings; Setting -x changes multiple parameters unless overridden. 
        pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref). 
        ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref).
        intractg: -B9 -O16 -L5  (intra-species contigs to ref).
      type:
      - 'null'
      - name: read_type
        type: enum
        symbols:
        - pacbio
        - ont2d
        - intractg
      inputBinding:
        prefix: -x
        position: 4
        shellQuote: false
      sbg:category: BWA Scoring options
    - id: reference_index_tar
      label: Reference Index TAR
      doc: Reference fasta file with its BWA index files packed in a TAR archive.
      type: File
      sbg:category: Input files
      sbg:fileTypes: TAR
    - id: mark_shorter
      label: Mark shorter
      doc: Mark shorter split hits as secondary.
      type: boolean?
      inputBinding:
        prefix: -M
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
    - id: speficy_distribution_parameters
      label: Specify distribution parameters
      doc: |-
        Specify the mean, standard deviation (10% of the mean if absent), max (4 sigma from the mean if absent), and min of the insert size distribution. 
        FR orientation only. 
        This array can have maximum of four values, where the first two should be specified as FLOAT and the last two as INT.
      type: float[]?
      inputBinding:
        prefix: -I
        position: 4
        valueFrom: |-
          ${
              var out = "";
              for (var i = 0; i < [].concat(self).length; i++ ){
                  out += " -I" + [].concat(self)[i];
              }    
              return out
          }
        separate: false
        itemSeparator: ' -I'
        shellQuote: false
      sbg:category: BWA Input/output options
    - id: minimum_output_score
      label: Minimum alignment score for a read to be output in SAM/BAM
      doc: Minimum alignment score for a read to be output in SAM/BAM.
      type: int?
      inputBinding:
        prefix: -T
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
      sbg:toolDefaultValue: '30'
    - id: output_format
      label: Output format
      doc: Coordinate sorted BAM file (option BAM) is the default output.
      type:
      - 'null'
      - name: output_format
        type: enum
        symbols:
        - SAM
        - BAM
        - CRAM
        - Queryname Sorted BAM
        - Queryname Sorted SAM
      sbg:category: Execution
      sbg:toolDefaultValue: Coordinate Sorted BAM
    - id: skip_mate_rescue
      label: Skip mate rescue
      doc: Skip mate rescue.
      type: boolean?
      inputBinding:
        prefix: -S
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
    - id: skip_seeds
      label: Skip seeds
      doc: Skip seeds with more than a given number (INT) of occurrences.
      type: int?
      inputBinding:
        prefix: -c
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '500'
    - id: output_name
      label: Output alignements file name
      doc: Name for the output alignments (SAM, BAM, or CRAM) file.
      type: string?
      sbg:category: Configuration
    - id: minimum_seed_length
      label: Minimum seed length
      doc: Minimum seed length for BWA MEM.
      type: int?
      inputBinding:
        prefix: -k
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '19'
    - id: gap_open_penalties
      label: Gap open penalties
      doc: |-
        Gap open penalties for deletions and insertions. 
        This array can't have more than two values.
      type: int[]?
      inputBinding:
        prefix: -O
        position: 4
        separate: false
        itemSeparator: ','
        shellQuote: false
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '[6,6]'
    - id: rg_median_fragment_length
      label: Median fragment length
      doc: Specify the median fragment length for RG line.
      type: string?
      sbg:category: BWA Read Group Options
    - id: mismatch_penalty
      label: Mismatch penalty
      doc: Penalty for a mismatch.
      type: int?
      inputBinding:
        prefix: -B
        position: 4
        shellQuote: false
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '4'
    - id: output_alignments
      label: Output alignments
      doc: Output all alignments for SE or unpaired PE.
      type: boolean?
      inputBinding:
        prefix: -a
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
    - id: discard_exact_matches
      label: Discard exact matches
      doc: Discard full-length exact matches.
      type: boolean?
      inputBinding:
        prefix: -e
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
    - id: rg_platform_unit_id
      label: Platform unit ID
      doc: |-
        Specify the platform unit (lane/slide) for RG line - An identifier for lanes (Illumina), or for slides (SOLiD) in the case that a library was split and ran over multiple lanes on the flow cell or slides.
      type: string?
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Inferred from metadata
    - id: mapQ_of_suplementary
      label: Don't modify mapQ
      doc: Don't modify mapQ of supplementary alignments.
      type: boolean?
      inputBinding:
        prefix: -q
        position: 4
        shellQuote: false
    - id: rg_sample_id
      label: Sample ID
      doc: |-
        Specify the sample ID for RG line - A human readable identifier for a sample or specimen, which could contain some metadata information. A sample or specimen is material taken from a biological entity for testing, diagnosis, propagation, treatment, or research purposes, including but not limited to tissues, body fluids, cells, organs, embryos, body excretory products, etc.
      type: string?
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Inferred from metadata
    - id: rg_data_submitting_center
      label: Data submitting center
      doc: Specify the data submitting center for RG line.
      type: string?
      sbg:category: BWA Read Group Options
    - id: discard_chain_length
      label: Discard chain length
      doc: Discard a chain if seeded bases are shorter than a given number (INT).
      type: int?
      inputBinding:
        prefix: -W
        position: 4
        shellQuote: false
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '0'
    - id: split_alignment_primary
      label: Split alignment - smallest coordinate as primary
      doc: for split alignment, take the alignment with the smallest coordinate as
        primary.
      type: boolean?
      inputBinding:
        prefix: '-5'
        position: 4
        shellQuote: false
    - id: append_comment
      label: Append comment
      doc: Append FASTA/FASTQ comment to the output file.
      type: boolean?
      inputBinding:
        prefix: -C
        position: 4
        shellQuote: false
      sbg:category: BWA Input/output options
    - id: read_group_header
      label: Read group header
      doc: |-
        Read group header line such as '@RG\tID:foo\tSM:bar'.  This value takes precedence over per-attribute parameters.
      type: string?
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Constructed from per-attribute parameters or inferred
        from metadata.
    - id: ignore_default_rg_id
      label: Ignore default RG ID
      doc: Ignore default RG ID ('1').
      type: boolean?
      sbg:category: BWA Read Group Options
    - id: fasta_index
      label: Fasta Index file for CRAM output
      doc: |-
        Fasta index file is required for CRAM output when no PCR Deduplication is selected.
      type: File?
      inputBinding:
        position: 4
        valueFrom: "${\n    return \"\";\n}"
        shellQuote: false
      sbg:category: Input files
      sbg:fileTypes: FAI

    outputs:
    - id: aligned_reads
      label: Aligned SAM/BAM
      doc: Aligned reads.
      type: File?
      secondaryFiles:
      - .bai
      - ^.bai
      - .crai
      - ^.crai
      outputBinding:
        glob: "${ \n    return [\"*.sam\", \"*.bam\", \"*.cram\"] \n}"
        outputEval: |-
          ${  
              /// Set metadata from input parameters, metadata or default value

              function flatten(files){
                  var a = []
                  for(var i=0;i<files.length;i++){
                      if(files[i]){
                          if(files[i].constructor == Array) a = a.concat(flatten(files[i]));
                          else a = a.concat(files[i]);}}
                  var b = a.filter(function (el) {return el != null});
                  return b;
              }
              function sharedStart(array){
                  var A= array.concat().sort(), 
                  a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                  while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                  return a1.substring(0, i);
              }
              /// Key-setting functions
              // Reference genome 
              var add_metadata_key_reference_genome = function(self, inputs) {
                  var reference_file = inputs.reference_index_tar.basename;
                  var ref_list = reference_file.split('.');
                  var  a = '';
                  a = ref_list.pop();
                  a = ref_list.pop();
                  a = ref_list.pop();
                  a = ref_list.pop(); // strip '.bwa-mem2-2.1-index-archive.tar'
                  return ref_list.join('.');
              };
              // Platform 
              var add_metadata_key_platform = function(self, inputs) {
                  /// Set platform from input parameters/input metadata/default value
                  var platform = '';
                  var pl = '';
                  // Find PL from header
                  if (inputs.read_group_header){
                      var header = inputs.read_group_header;
                      header = header.split("'").join("") //remove single quotes
                      var a = header.split('\\t');
                      for (var i = 0; i < a.length; i++){ //find PL field
                          if (a[i].includes("PL:")) pl= a[i];
                          else;
                      }}
                  else;
                  
                  if (pl) platform = pl.split(':')[1];
                  else if (inputs.rg_platform) platform = inputs.rg_platform;
                  else if (read_metadata.platform) platform = read_metadata.platform;
                  else platform = 'Illumina';
                  
                  return platform
              };
              // Sample ID 
              var add_metadata_key_sample_id = function(self, inputs) {
                  /// Set sample ID from input parameters/input metadata/default value from input reads file names
                  var sample_id = '';
                  var sm = '';
                  // Find SM from header
                  if (inputs.read_group_header){
                      var header = inputs.read_group_header;
                      header = header.split("'").join("") //remove single quotes
                      var a = header.split('\\t');
                      for (var i = 0; i < a.length; i++){ //find SM field
                          if (a[i].includes("SM:")) var sm= a[i];
                          else;
                      }}
                  else;
                  
                  if (sm) sample_id = sm.split(':')[1];
                  else if (inputs.rg_sample_id) sample_id = inputs.rg_sample_id;
                  else if (read_metadata.sample_id) sample_id = read_metadata.sample_id;
                  else {
                      var read_names = [];
                      var files1 = [].concat(inputs.input_reads);
                      var files=flatten(files1);
                      
                      for (var i=0;i<files.length;i++) {
                          var file_ext=files[i].nameext;
                          var file_base=files[i].basename;
                          
                          if (file_ext === '.gz' || file_ext === '.GZ')
                              file_base = file_base.slice(0, -3);
                              file_ext= '.'+ file_base.split('.').pop();
                          if (file_ext === '.fq' || file_ext === '.FQ')
                              file_base = file_base.slice(0, -3);
                          if (file_ext === '.fastq' || file_ext === '.FASTQ')
                              file_base = file_base.slice(0, -6);
                          
                          read_names.push(file_base.replace(/pe1|pe2|pe\.1|pe\.2|pe\_1|pe\_2|\_pe1|\_pe2|\_pe\.1|\_pe\.2|\_pe\_1|\_pe\_2|\.pe1|\.pe2|\.pe\.1|\.pe\.2|\.pe\_1|\.pe\_2/,''));
                        }
                        ////strip out any trailing dashes/dots/underscores...
                        var unique_prefix = sharedStart(read_names).replace( /\-$|\_$|\.$/, '');
                        var tmp_prefix = unique_prefix.replace( /^\_|\.pe$|\.R$|\_pe$|\_R$/,'');
                        var final_prefix = tmp_prefix.replace( /^_\d(\d)?_/, '' );
                        
                        var fname=final_prefix;
                      sample_id = fname;
                  }
                  return sample_id
              };
              
             
              var files1 = [].concat(inputs.input_reads);
              var files=flatten(files1);
              var read_metadata = files[0].metadata;
              if (!read_metadata) read_metadata = [];
              
              self = inheritMetadata(self, files);

              for (var i = 0; i < self.length; i++) {
                  var out_metadata = {
                      'reference_genome': add_metadata_key_reference_genome(self[i], inputs),
                      'platform': add_metadata_key_platform(self[i], inputs),
                      'sample_id': add_metadata_key_sample_id(self[i], inputs)
                  };
                  self[i] = setMetadata(self[i], out_metadata);
              }

              return self;

          }
      sbg:fileTypes: SAM, BAM, CRAM
    - id: dups_metrics
      label: Sormadup metrics
      doc: Metrics file for biobambam mark duplicates
      type: File?
      outputBinding:
        glob: '*.sormadup_metrics.log'
      sbg:fileTypes: LOG

    baseCommand: []
    arguments:
    - prefix: ''
      position: -1
      valueFrom: |-
        ${
            /// Check number of input FASTQ files ///
            
            function flatten(files){
            var a = []
            for(var i=0;i<files.length;i++){
                if(files[i]){
                    if(files[i].constructor == Array) a = a.concat(flatten(files[i]));
                    else a = a.concat(files[i])}}
                var b = a.filter(function (el) {return el != null})
                return b
            }
            
            var files1 = [].concat(inputs.input_reads);
            var in_reads=flatten(files1);
            
            if ( in_reads.length > 2 ) return 'ERROR: Number of input FASTQ files needs to be one (if single-end/interleaved file) or two (if paired-end files)';
            else return '';
        }
      shellQuote: false
    - prefix: ''
      position: 0
      valueFrom: |-
        ${
            var cmd = "/bin/bash -c \"";
            return cmd + " export REF_CACHE=${PWD} && ";
        }
      shellQuote: false
    - prefix: ''
      position: 1
      valueFrom: |-
        ${
            /// Unpack Reference TAR archive ///
            
            var in_index=[].concat(inputs.reference_index_tar)[0];
            var reference_file = in_index.basename;
            return 'tar -tvf ' + reference_file + ' 1>&2 && tar -xf ' + reference_file + ' && ';
            
        }
      shellQuote: false
    - prefix: ''
      position: 2
      valueFrom: bwa mem
      shellQuote: false
    - prefix: ''
      position: 5
      valueFrom: |-
        ${
            /// Set RG header ///

            function add_param(key, val) {
                if (!val) return;
                param_list.push(key + ':' + val);}
                
            function flatten(files){
                var a = [];
                for(var i=0;i<files.length;i++){
                    if(files[i]){
                        if(files[i].constructor == Array) a = a.concat(flatten(files[i]));
                        else a = a.concat(files[i]);}}
                var b = a.filter(function (el) {return el != null;});
                return b;}
                
            function sharedStart(array){
                var A= array.concat().sort(), 
                a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                return a1.substring(0, i);}

            
            /// If it exists - return input read group header from input parameter
            if (inputs.read_group_header) return '-R ' + inputs.read_group_header;

            // Flatten input reads
            var in_reads1 = [].concat(inputs.input_reads);
            var in_reads = flatten(in_reads1)
            var input_1=in_reads[0];

            var param_list = [];
            //Read metadata for input reads
            var read_metadata = input_1.metadata;
            if (!read_metadata) read_metadata = [];

            // Set CN
            if (inputs.rg_data_submitting_center) add_param('CN', inputs.rg_data_submitting_center);
            else if ('data_submitting_center' in read_metadata) add_param('CN', read_metadata.data_submitting_center);
            else;

            // Set LB
            if (inputs.rg_library_id) add_param('LB', inputs.rg_library_id);
            else if ('library_id' in read_metadata) add_param('LB', read_metadata.library_id);
            else;

            // Set PI
            if (inputs.rg_median_fragment_length) add_param('PI', inputs.rg_median_fragment_length);
            else;

            // Set PL (default Illumina)
            var rg_platform = '';
            if (inputs.rg_platform) add_param('PL', inputs.rg_platform);
            else if ('platform' in read_metadata) {
                if (read_metadata.platform == 'HiSeq X Ten') rg_platform = 'Illumina';
                else rg_platform = read_metadata.platform;
                add_param('PL', rg_platform);}
            else add_param('PL', 'Illumina');

            // Set PU
            if (inputs.rg_platform_unit_id) add_param('PU', inputs.rg_platform_unit_id);
            else if ('platform_unit_id' in read_metadata) add_param('PU', read_metadata.platform_unit_id);
            else;
            
            // Set RG_ID
            var folder = input_1.path.split('/').slice(-2,-1).toString();
            var suffix = "_s";
            
            if (inputs.rg_id) add_param('ID', inputs.rg_id);
            else if (folder.indexOf(suffix, folder.length - suffix.length) !== -1){/// Set unique RG_ID when in scatter mode
                var rg = folder.split("_").slice(-2)[0];
                if (parseInt(rg)) add_param('ID', rg);
                else add_param('ID', 1);}
            else  add_param('ID', 1);

            // Set SM from input/metadata/filename
            if (inputs.rg_sample_id) add_param('SM', inputs.rg_sample_id);
            else if ('sample_id' in read_metadata) add_param('SM', read_metadata.sample_id);
            else {
                var read_names = [];
                for (var i=0;i<in_reads.length;i++) {
                    var file_ext=in_reads[i].nameext;
                    var file_base=in_reads[i].basename;
                    
                    if (file_ext === '.gz' || file_ext === '.GZ')
                        file_base = file_base.slice(0, -3);
                        file_ext= '.'+ file_base.split('.').pop();
                    if (file_ext === '.fq' || file_ext === '.FQ')
                        file_base = file_base.slice(0, -3);
                    if (file_ext === '.fastq' || file_ext === '.FASTQ')
                        file_base = file_base.slice(0, -6);
                    
                    read_names.push(file_base.replace(/pe1|pe2|pe\.1|pe\.2|pe\_1|pe\_2|\_pe1|\_pe2|\_pe\.1|\_pe\.2|\_pe\_1|\_pe\_2|\.pe1|\.pe2|\.pe\.1|\.pe\.2|\.pe\_1|\.pe\_2/,''));}
                  
                ////strip out any trailing dashes/dots/underscores...
                var unique_prefix = sharedStart(read_names).replace( /\-$|\_$|\.$/, '');
                var tmp_prefix = unique_prefix.replace( /^\_|\.pe$|\.R$|\_pe$|\_R$/,'');
                var final_prefix = tmp_prefix.replace( /^_\d(\d)?_/, '' );
              
                var sample_id=final_prefix;
                add_param('SM', sample_id);
            };
            
            if (!inputs.ignore_default_rg_id) {
              return "-R '@RG\\t" + param_list.join('\\t') + "'";
            } else {
              return '';
            }

        }
      shellQuote: false
    - prefix: -t
      position: 6
      valueFrom: |-
        ${
            /// Set BWA2 threads ///

            var  MAX_THREADS = 36;
            var  suggested_threads = 8;
            var threads  = 0;
          
            if (inputs.threads) threads = inputs.threads;
            else if (inputs.wgs_hg38_mode_threads) {
                var ref_name = inputs.reference_index_tar.basename;
                if (ref_name.search('38') >= 0) threads = inputs.wgs_hg38_mode_threads;
                else threads = MAX_THREADS;
            } else threads = suggested_threads;
            
            return threads;
        }
      shellQuote: false
    - prefix: ''
      position: 14
      valueFrom: |-
        ${
            /// Extract common prefix for Index files ///
            
            var reference_tar = [].concat(inputs.reference_index_tar)[0];
            
            var prefix = "$(tar -tf " + reference_tar.basename + " --wildcards '*.bwt' | rev | cut -c 5- | rev)";
            return prefix;

        }
      shellQuote: false
    - prefix: ''
      position: 116
      valueFrom: |-
        ${
            ///  BIOBAMBAM2  ///
              
             // Get shared start and flatten input reads
            function sharedStart(array){
                var A= array.concat().sort(), 
                a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;
                while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;
                return a1.substring(0, i);
            }
            function flatten(files){
                var a = [];
                for(var i=0;i<files.length;i++){
                    if(files[i]){
                        if(files[i].constructor == Array) a = a.concat(flatten(files[i]));
                        else a = a.concat(files[i]);}}
                var b = a.filter(function (el) {return el != null;});
                return b;}
           
            var input_reads = [].concat(inputs.input_reads);
            var files=flatten(input_reads);

            // Set output file name
            var fname = '';
            
            /// from given prefix
            if (inputs.output_name) fname = inputs.output_name;
            /// from sample_id metadata
            else if (files[0].metadata && files[0].metadata['sample_id']) fname=files[0].metadata['sample_id'];
            /// from common prefix, and strip out any unnecessary characters
            else {
                var read_names = [];
                for (var i=0;i<files.length;i++) {
                    var file_ext=files[i].nameext;
                    var file_base=files[i].basename;
                    
                    if (file_ext === '.gz' || file_ext === '.GZ')
                        file_base = file_base.slice(0, -3);
                        file_ext= '.'+ file_base.split('.').pop();
                    if (file_ext === '.fq' || file_ext === '.FQ')
                        file_base = file_base.slice(0, -3);
                    if (file_ext === '.fastq' || file_ext === '.FASTQ')
                        file_base = file_base.slice(0, -6);
                    
                    read_names.push(file_base.replace(/pe1|pe2|pe\.1|pe\.2|pe\_1|pe\_2|\_pe1|\_pe2|\_pe\.1|\_pe\.2|\_pe\_1|\_pe\_2|\.pe1|\.pe2|\.pe\.1|\.pe\.2|\.pe\_1|\.pe\_2/,''));
                      
                  }
                  ////strip out any trailing dashes/dots/underscores...
                  var unique_prefix = sharedStart(read_names).replace( /\-$|\_$|\.$/, '');
                  var tmp_prefix = unique_prefix.replace( /^\_|\.pe$|\.R$|\_pe$|\_R$/,'');
                  var final_prefix = tmp_prefix.replace( /^_\d(\d)?_/, '' );
                  
                  fname=final_prefix;}


            // Read number of threads if defined
            var threads = 0;
            var MAX_THREADS = 0;
            var ref_name = '';
            if (inputs.threads) threads = inputs.threads;
            else if (inputs.wgs_hg38_mode_threads) {
                MAX_THREADS = 36;
                ref_name = inputs.reference_index_tar.basename;
                if (ref_name.search('38') >= 0) threads = inputs.wgs_hg38_mode_threads;
                else threads = MAX_THREADS;
                } 
            else threads = 8;

            var tool = '';
            var dedup = '';
            if (inputs.deduplication == "MarkDuplicates") {
                tool = 'bamsormadup';
                dedup = ' markduplicates=1';
            } else {
                if (inputs.output_format == 'CRAM') tool = 'bamsort index=0';
                else tool = 'bamsort index=1';
                if (inputs.deduplication == "RemoveDuplicates") dedup = ' rmdup=1';
                else dedup = '';
            }
            var sort_path = tool + dedup;

            var indexfilename = '';
            var out_format = '';
            var extension  = '';
            // Coordinate Sorted BAM is default
            if (inputs.output_format == 'CRAM') {
                out_format = ' outputformat=cram SO=coordinate';
                ref_name = inputs.reference_index_tar.basename.split('.tar')[0];
                out_format += ' reference=' + ref_name;
                if (sort_path != 'bamsort index=0') indexfilename = ' indexfilename=' + fname + '.cram.crai';
                extension = '.cram';
            } else if (inputs.output_format == 'SAM') {
                out_format = ' outputformat=sam SO=coordinate';
                extension = '.sam';
            } else if (inputs.output_format == 'Queryname Sorted BAM') {
                out_format = ' outputformat=bam SO=queryname';
                extension = '.bam';
            } else if (inputs.output_format == 'Queryname Sorted SAM') {
                out_format = ' outputformat=sam SO=queryname';
                extension = '.sam';
            } else {
                out_format = ' outputformat=bam SO=coordinate';
                indexfilename = ' indexfilename=' + fname + '.bam.bai';
                extension = '.bam';
            }
            var cmd = " | " + sort_path + " threads=" + threads + " level=1 tmplevel=-1 inputformat=sam";
            cmd += out_format;
            cmd += indexfilename;
            // capture metrics file
            cmd += " M=" + fname + ".sormadup_metrics.log";

            if (inputs.output_format == 'SAM') cmd = '';
            
            return cmd + ' > ' + fname + extension;
            
        }
      separate: false
      shellQuote: false
    - prefix: ''
      position: 10004
      valueFrom: |-
        ${
            /// Get pipe status ///
            
            var  cmd = ";declare -i pipe_statuses=(\\${PIPESTATUS[*]});len=\\${#pipe_statuses[@]};declare -i tot=0;echo \\${pipe_statuses[*]};for (( i=0; i<\\${len}; i++ ));do if [ \\${pipe_statuses[\\$i]} -ne 0 ];then tot=\\${pipe_statuses[\\$i]}; fi;done;if [ \\$tot -ne 0 ]; then >&2 echo Error in piping. Pipe statuses: \\${pipe_statuses[*]};fi; if [ \\$tot -ne 0 ]; then false;fi\"";
            return cmd;
        }
      shellQuote: false
    id: nens/bwa-0-7-15-cwl1-0-demo/bwa-mem-bundle-0-7-15/21
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Genomics
    - Alignment
    - CWL1.0
    sbg:cmdPreview: |-
      /bin/bash -c " export REF_CACHE=${PWD} ;  tar -tvf reference.HG38.fasta.gz.tar 1>&2; tar -xf reference.HG38.fasta.gz.tar ;  bwa mem  -R '@RG\tID:1\tPL:Illumina\tSM:dnk_sample' -t 10  reference.HG38.fasta.gz  /path/to/LP6005524-DNA_C01_lane_7.sorted.converted.filtered.pe_2.gz /path/to/LP6005524-DNA_C01_lane_7.sorted.converted.filtered.pe_1.gz  | bamsormadup threads=8 level=1 tmplevel=-1 inputformat=sam outputformat=cram SO=coordinate reference=reference.HG38.fasta.gz indexfilename=LP6005524-DNA_C01_lane_7.sorted.converted.filtered.cram.crai M=LP6005524-DNA_C01_lane_7.sorted.converted.filtered.sormadup_metrics.log > LP6005524-DNA_C01_lane_7.sorted.converted.filtered.cram  ;declare -i pipe_statuses=(\${PIPESTATUS[*]});len=\${#pipe_statuses[@]};declare -i tot=0;echo \${pipe_statuses[*]};for (( i=0; i<\${len}; i++ ));do if [ \${pipe_statuses[\$i]} -ne 0 ];then tot=\${pipe_statuses[\$i]}; fi;done;if [ \$tot -ne 0 ]; then >&2 echo Error in piping. Pipe statuses: \${pipe_statuses[*]};fi; if [ \$tot -ne 0 ]; then false;fi"
    sbg:content_hash: a4965586211232dc4651281d3de154eac59adbbe47becb0c3a5f73560b751f560
    sbg:contributors:
    - nens
    - ana_stankovic
    - uros_sipetic
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1555689212
    sbg:expand_workflow: false
    sbg:id: h-167b8029/h-3f6bacf5/h-d72ab5d5/0
    sbg:image_url:
    sbg:latestRevision: 21
    sbg:license: |-
      BWA: GNU Affero General Public License v3.0, MIT License; Biobambam2: GNU General Public License v3.0
    sbg:links:
    - id: http://bio-bwa.sourceforge.net/
      label: Homepage
    - id: https://github.com/lh3/bwa
      label: Source code
    - id: http://bio-bwa.sourceforge.net/bwa.shtml
      label: Wiki
    - id: http://sourceforge.net/projects/bio-bwa/
      label: Download
    - id: http://arxiv.org/abs/1303.3997
      label: Publication
    - id: http://www.ncbi.nlm.nih.gov/pubmed/19451168
      label: Publication BWA Algorithm
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1611175341
    sbg:project: nens/bwa-0-7-15-cwl1-0-demo
    sbg:projectName: BWA 0.7.15 CWL1.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 21
    sbg:revisionNotes: added ignore_rg_id
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1555689212
      sbg:revision: 0
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/1
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1556035789
      sbg:revision: 1
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/3
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1556037315
      sbg:revision: 2
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/4
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1556192655
      sbg:revision: 3
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/5
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1556193727
      sbg:revision: 4
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/6
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000453
      sbg:revision: 5
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/9
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558002186
      sbg:revision: 6
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/10
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558021975
      sbg:revision: 7
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/12
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558023132
      sbg:revision: 8
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/13
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558085159
      sbg:revision: 9
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/15
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558349205
      sbg:revision: 10
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/16
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351490
      sbg:revision: 11
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/17
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558427784
      sbg:revision: 12
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/18
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558441939
      sbg:revision: 13
      sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/22
    - sbg:modifiedBy: ana_stankovic
      sbg:modifiedOn: 1579532841
      sbg:revision: 14
      sbg:revisionNotes: Bug fix for CRAM output with no PCR deduplication
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1581075318
      sbg:revision: 15
      sbg:revisionNotes: dev - v25; var added
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1581350490
      sbg:revision: 16
      sbg:revisionNotes: |-
        Add platform read group to the BAM even when no_rg_information parameter is specified, based on the input BAM platform metadata.
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1581359515
      sbg:revision: 17
      sbg:revisionNotes: Remove the default PL RG bit
    - sbg:modifiedBy: ana_stankovic
      sbg:modifiedOn: 1592998681
      sbg:revision: 18
      sbg:revisionNotes: Updated JS to assign a unique Read group ID when the tool
        is scattered
    - sbg:modifiedBy: ana_stankovic
      sbg:modifiedOn: 1609141711
      sbg:revision: 19
      sbg:revisionNotes: |-
        JavaScript cleanup; Default setting of Platform and Sample ID; Description update
    - sbg:modifiedBy: ana_stankovic
      sbg:modifiedOn: 1609169898
      sbg:revision: 20
      sbg:revisionNotes: filter_out_secondary_alignments parameter removed
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1611175341
      sbg:revision: 21
      sbg:revisionNotes: added ignore_rg_id
    sbg:sbgMaintained: false
    sbg:toolAuthor: Heng Li
    sbg:toolkit: BWA
    sbg:toolkitVersion: 0.7.15
    sbg:validationErrors: []
  out:
  - id: aligned_reads
  - id: dups_metrics
  sbg:x: -334.3309020996094
  sbg:y: 257.5992736816406
- id: gatk_mergebamalignment_4_1_0_0
  label: GATK MergeBamAlignment
  in:
  - id: add_mate_cigar
    default: 'true'
  - id: in_alignments
    valueFrom: $([self])
    source:
    - samtools_view_1_9_cwl1_0/out_alignments
  - id: aligner_proper_pair_flags
    default: true
  - id: attributes_to_retain
    default:
    - X0
  - id: clip_adapters
    default: 'false'
  - id: expected_orientations
    default:
    - FR
  - id: max_insertions_or_deletions
    default: -1
  - id: max_records_in_ram
    default: 2000000
  - id: paired_run
    default: 'true'
  - id: primary_alignment_strategy
    default: MostDistant
  - id: program_group_command_line
    default: '"bwa mem -K 100000000 -p -v 3 -t 16 -Y ref_fasta"'
  - id: program_group_name
    default: bwamem
  - id: program_group_version
    default: 0.7.15
  - id: program_record_id
    default: bwamem
  - id: in_reference
    source: in_reference
  - id: sort_order
    default: unsorted
  - id: unmap_contaminant_reads
    default: true
  - id: unmapped_bam
    source: in_alignments
  - id: unmapped_read_strategy
    default: COPY_TO_TAG
  - id: validation_stringency
    default: SILENT
  scatter:
  - in_alignments
  - unmapped_bam
  scatterMethod: dotproduct
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK MergeBamAlignment
    doc: |-
      The **GATK MergeBamAlignment** tool is used for merging BAM/SAM alignment info from a third-party aligner with the data in an unmapped BAM file, producing a third BAM file that has alignment data (from the aligner) and all the remaining data from the unmapped BAM.

      Many alignment tools still require FASTQ format input. The unmapped BAM may contain useful information that will be lost in the conversion to FASTQ (meta-data like sample alias, library, barcodes, etc... as well as read-level tags.) This tool takes an unaligned BAM with meta-data, and the aligned BAM produced by calling [SamToFastq](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_SamToFastq.php) and then passing the result to an aligner. It produces a new SAM file that includes all aligned and unaligned reads and also carries forward additional read attributes from the unmapped BAM (attributes that are otherwise lost in the process of converting to FASTQ). The resulting file will be valid for use by Picard and GATK tools. The output may be coordinate-sorted, in which case the tags, NM, MD, and UQ will be calculated and populated, or query-name sorted, in which case the tags will not be calculated or populated [1].

      *A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*

      ###Common Use Cases

      * The **GATK MergeBamAlignment** tool requires a SAM or BAM file on its **Aligned BAM/SAM file** (`--ALIGNED_BAM`) input, original SAM or BAM file of unmapped reads, which must be in queryname order on its **Unmapped BAM/SAM file** (`--UNMAPPED_BAM`) input and a reference sequence on its **Reference** (`--REFERENCE_SEQUENCE`) input. The tool generates a single BAM/SAM file on its **Output merged BAM/SAM file** output.

      * Usage example:

      ```
      gatk MergeBamAlignment \\
            --ALIGNED_BAM aligned.bam \\
            --UNMAPPED_BAM unmapped.bam \\
            --OUTPUT merged.bam \\
            --REFERENCE_SEQUENCE reference_sequence.fasta
      ```

      ###Changes Introduced by Seven Bridges

      * The output file name will be prefixed using the **Output prefix** parameter. In case **Output prefix** is not provided, output prefix will be the same as the Sample ID metadata from **Input SAM/BAM file**, if the Sample ID metadata exists. Otherwise, output prefix will be inferred from the **Input SAM/BAM file** filename. This way, having identical names of the output files between runs is avoided. Moreover,  **merged** will be added before the extension of the output file name. 

      * The user has a possibility to specify the output file format using the **Output file format** argument. Otherwise, the output file format will be the same as the format of the input aligned file.

      ###Common Issues and Important Notes

      * Note:  This is not a tool for taking multiple BAM/SAM files and creating a bigger file by merging them. For that use-case, see [MergeSamFiles](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_MergeSamFiles.php).

      ###Performance Benchmarking

      Below is a table describing runtimes and task costs of **GATK MergeBamAlignment** for a couple of different samples, executed on the AWS cloud instances:

      | Experiment type |  Aligned BAM/SAM size |  Unmapped BAM/SAM size | Duration |  Cost | Instance (AWS) | 
      |:--------------:|:------------:|:--------:|:-------:|:---------:|:----------:|:------:|:------:|------:|
      |     RNA-Seq     |  1.4 GB |  1.9 GB |   9min   | ~0.06$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     |  4.0 GB |  5.7 GB |   20min   | ~0.13$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 6.6 GB | 9.5 GB |  32min  | ~0.21$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 13 GB | 19 GB |  1h 4min  | ~0.42$ | c4.2xlarge (8 CPUs) |

      *Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*

      ###References

      [1] [GATK MergeBamAlignment](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_MergeBamAlignment.php)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: "${\n    return inputs.cpu_per_job ? inputs.cpu_per_job : 1;\n}"
      ramMin: |-
        ${
            var memory = 4096;
            if (inputs.memory_per_job) 
            {
                memory = inputs.memory_per_job;
            }
            if (inputs.memory_overhead_per_job)
            {
                memory += inputs.memory_overhead_per_job;
            }
            return memory;
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

    inputs:
    - id: add_mate_cigar
      label: Add mate CIGAR
      doc: Adds the mate CIGAR tag (MC) if true, does not if false.
      type:
      - 'null'
      - name: add_mate_cigar
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --ADD_MATE_CIGAR
        position: 4
        shellQuote: false
      sbg:altPrefix: -MC
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: add_pg_tag_to_reads
      label: Add PG tag to reads
      doc: Add PG tag to each read in a SAM or BAM.
      type:
      - 'null'
      - name: add_pg_tag_to_reads
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --ADD_PG_TAG_TO_READS
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: in_alignments
      label: Aligned BAM/SAM file
      doc: |-
        SAM or BAM file(s) with alignment data. Cannot be used in conjuction with argument(s) READ1_ALIGNED_BAM (R1_ALIGNED) READ2_ALIGNED_BAM (R2_ALIGNED).
      type: File[]
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              var arr = [].concat(inputs.in_alignments);
              if (arr.length == 1) 
              {
                  return "--ALIGNED_BAM " + arr[0].path;
              }
              else
              {
                  var pe_1 = [];
                  var pe_2 = [];
                  var se = [];
                  for (var i in arr)
                  {
                      if (arr[i].metadata && arr[i].metadata.paired_end && arr[i].metadata.paired_end == 1)
                      {
                          pe_1.push(arr[i].path);
                      }
                      else if (arr[i].metadata && arr[i].metadata.paired_end && arr[i].metadata.paired_end == 2)
                      {
                          pe_2.push(arr[i].path);
                      }
                      else
                      {
                          se.push(arr[i].path);
                      }
                  }
                  
                  if (se.length > 0) 
                  {
                      return "--ALIGNED_BAM " + se.join(" --ALIGNED_BAM ");
                  } 
                  else if (pe_1.length > 0 && pe_2.length > 0 && pe_1.length == pe_2.length) 
                  {
                      return "--READ1_ALIGNED_BAM " + pe_1.join(' --READ1_ALIGNED_BAM ') + " --READ2_ALIGNED_BAM " + pe_2.join(' --READ2_ALIGNED_BAM ');
                  } 
                  else 
                  {
                      return "";
                  }
                      
              }
          }
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:fileTypes: BAM, SAM
      sbg:toolDefaultValue: 'null'
    - id: aligned_reads_only
      label: Aligned reads only
      doc: Whether to output only aligned reads.
      type: boolean?
      inputBinding:
        prefix: --ALIGNED_READS_ONLY
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: aligner_proper_pair_flags
      label: Aligner proper pair flags
      doc: |-
        Use the aligner's idea of what a proper pair is rather than computing in this program.
      type: boolean?
      inputBinding:
        prefix: --ALIGNER_PROPER_PAIR_FLAGS
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: attributes_to_remove
      label: Attributes to remove
      doc: |-
        Attributes from the alignment record that should be removed when merging. This overrides ATTRIBUTES_TO_RETAIN if they share common tags.
      type: string[]?
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--ATTRIBUTES_TO_REMOVE', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: attributes_to_retain
      label: Attributes to retain
      doc: |-
        Reserved alignment attributes (tags starting with X, Y, or Z) that should be brought over from the alignment data when merging.
      type: string[]?
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--ATTRIBUTES_TO_RETAIN', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: attributes_to_reverse
      label: Attributes to reverse
      doc: Attributes on negative strand reads that need to be reversed.
      type: string[]?
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--ATTRIBUTES_TO_REVERSE', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:altPrefix: -RV
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '[OQ,U2]'
    - id: attributes_to_reverse_complement
      label: Attributes to reverse complement
      doc: Attributes on negative strand reads that need to be reverse complemented.
      type: string[]?
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--ATTRIBUTES_TO_REVERSE_COMPLEMENT', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:altPrefix: -RC
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '[E2,SQ]'
    - id: clip_adapters
      label: Clip adapters
      doc: Whether to clip adapters where identified.
      type:
      - 'null'
      - name: clip_adapters
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --CLIP_ADAPTERS
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: clip_overlapping_reads
      label: Clip overlapping reads
      doc: |-
        For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate.
      type:
      - 'null'
      - name: clip_overlapping_reads
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --CLIP_OVERLAPPING_READS
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: compression_level
      label: Compression level
      doc: Compression level for all compressed files created (e.g. BAM and VCF).
      type: int?
      inputBinding:
        prefix: --COMPRESSION_LEVEL
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
    - id: create_index
      label: Create index
      doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
      type: boolean?
      inputBinding:
        prefix: --CREATE_INDEX
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: expected_orientations
      label: Expected orientations
      doc: |-
        The expected orientation of proper read pairs. Replaces JUMP_SIZE. Cannot be used in conjuction with argument(s) JUMP_SIZE (JUMP).
      type: string[]?
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--EXPECTED_ORIENTATIONS', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:altPrefix: -ORIENTATIONS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: include_secondary_alignments
      label: Include secondary alignments
      doc: If false, do not write secondary alignments to output.
      type:
      - 'null'
      - name: include_secondary_alignments
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --INCLUDE_SECONDARY_ALIGNMENTS
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: is_bisulfite_sequence
      label: Is bisulfite sequence
      doc: Whether the lane is bisulfite sequence (used when calculating the NM tag).
      type: boolean?
      inputBinding:
        prefix: --IS_BISULFITE_SEQUENCE
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: jump_size
      label: Jump size
      doc: |-
        The expected jump size (required if this is a jumping library). Deprecated. Use EXPECTED_ORIENTATIONS instead. Cannot be used in conjuction with argument(s) EXPECTED_ORIENTATIONS (ORIENTATIONS).
      type: int?
      inputBinding:
        prefix: --JUMP_SIZE
        position: 4
        shellQuote: false
      sbg:altPrefix: -JUMP
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: matching_dictionary_tags
      label: Matching dictionary tags
      doc: |-
        List of Sequence Records tags that must be equal (if present) in the reference dictionary and in the aligned file. Mismatching tags will cause an error if in this list, and a warning otherwise.
      type: string[]?
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--MATCHING_DICTIONARY_TAGS', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '[M5,LN]'
    - id: max_insertions_or_deletions
      label: Max insertions or deletions
      doc: |-
        The maximum number of insertions or deletions permitted for an alignment to be included. Alignments with more than this many insertions or deletions will be ignored. Set to -1 to allow any number of insertions or deletions.
      type: int?
      inputBinding:
        prefix: --MAX_INSERTIONS_OR_DELETIONS
        position: 4
        shellQuote: false
      sbg:altPrefix: -MAX_GAPS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '1'
    - id: max_records_in_ram
      label: Max records in RAM
      doc: |-
        When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
      type: int?
      inputBinding:
        prefix: --MAX_RECORDS_IN_RAM
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
    - id: memory_overhead_per_job
      label: Memory overhead per job
      doc: |-
        This input allows a user to set the desired overhead memory when running a tool or adding it to a workflow. This amount will be added to the Memory per job in the Memory requirements section but it will not be added to the -Xmx parameter leaving some memory not occupied which can be used as stack memory (-Xmx parameter defines heap memory). This input should be defined in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      type: int?
      sbg:category: Platform Options
    - id: memory_per_job
      label: Memory per job
      doc: |-
        This input allows a user to set the desired memory requirement when running a tool or adding it to a workflow. This value should be propagated to the -Xmx parameter too.This input should be defined in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      type: int?
      sbg:category: Platform Options
    - id: min_unclipped_bases
      label: Min unclipped bases
      doc: |-
        If UNMAP_CONTAMINANT_READS is set, require this many unclipped bases or else the read will be marked as contaminant.
      type: int?
      inputBinding:
        prefix: --MIN_UNCLIPPED_BASES
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '32'
    - id: paired_run
      label: Paired run
      doc: DEPRECATED. This argument is ignored and will be removed.
      type:
      - 'null'
      - name: paired_run
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --PAIRED_RUN
        position: 4
        shellQuote: false
      sbg:altPrefix: -PE
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: primary_alignment_strategy
      label: Primary alignment strategy
      doc: |-
        Strategy for selecting primary alignment when the aligner has provided more than one alignment for a pair or fragment, and none are marked as primary, more than one is marked as primary, or the primary alignment is filtered out for some reason. For all strategies, ties are resolved arbitrarily. Possible values: { BestMapq (expects that multiple alignments will be correlated with HI tag, and prefers the pair of alignments with the largest MAPQ, in the absence of a primary selected by the aligner.) EarliestFragment (prefers the alignment which maps the earliest base in the read. Note that EarliestFragment may not be used for paired reads.) BestEndMapq (appropriate for cases in which the aligner is not pair-aware, and does not output the HI tag. It simply picks the alignment for each end with the highest MAPQ, and makes those alignments primary, regardless of whether the two alignments make sense together.) MostDistant (appropriate for a non-pair-aware aligner. Picks the alignment pair with the largest insert size. If all alignments would be chimeric, it picks the alignments for each end with the best MAPQ. ) }.
      type:
      - 'null'
      - name: primary_alignment_strategy
        type: enum
        symbols:
        - BestMapq
        - EarliestFragment
        - BestEndMapq
        - MostDistant
      inputBinding:
        prefix: --PRIMARY_ALIGNMENT_STRATEGY
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: BestMapq
    - id: program_group_command_line
      label: Program group command line
      doc: The command line of the program group (if not supplied by the aligned file).
      type: string?
      inputBinding:
        prefix: --PROGRAM_GROUP_COMMAND_LINE
        position: 4
        shellQuote: false
      sbg:altPrefix: -PG_COMMAND
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: program_group_name
      label: Program group name
      doc: The name of the program group (if not supplied by the aligned file).
      type: string?
      inputBinding:
        prefix: --PROGRAM_GROUP_NAME
        position: 4
        shellQuote: false
      sbg:altPrefix: -PG_NAME
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: program_group_version
      label: Program group version
      doc: The version of the program group (if not supplied by the aligned file).
      type: string?
      inputBinding:
        prefix: --PROGRAM_GROUP_VERSION
        position: 4
        shellQuote: false
      sbg:altPrefix: -PG_VERSION
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: program_record_id
      label: Program record id
      doc: The program group ID of the aligner (if not supplied by the aligned file).
      type: string?
      inputBinding:
        prefix: --PROGRAM_RECORD_ID
        position: 4
        shellQuote: false
      sbg:altPrefix: -PG
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: read1_trim
      label: Read1 trim
      doc: The number of bases trimmed from the beginning of read 1 prior to alignment.
      type: int?
      inputBinding:
        prefix: --READ1_TRIM
        position: 4
        shellQuote: false
      sbg:altPrefix: -R1_TRIM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: read2_trim
      label: Read2 trim
      doc: The number of bases trimmed from the beginning of read 2 prior to alignment.
      type: int?
      inputBinding:
        prefix: --READ2_TRIM
        position: 4
        shellQuote: false
      sbg:altPrefix: -R2_TRIM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: in_reference
      label: Reference
      doc: Reference sequence file.
      type: File
      secondaryFiles:
      - .fai
      - ^.dict
      inputBinding:
        prefix: --REFERENCE_SEQUENCE
        position: 4
        shellQuote: false
      sbg:altPrefix: -R
      sbg:category: Required Arguments
      sbg:fileTypes: FASTA, FA
    - id: sort_order
      label: Sort order
      doc: The order in which the merged reads should be output.
      type:
      - 'null'
      - name: sort_order
        type: enum
        symbols:
        - unsorted
        - queryname
        - coordinate
        - duplicate
        - unknown
      inputBinding:
        prefix: --SORT_ORDER
        position: 4
        shellQuote: false
      sbg:altPrefix: -SO
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: coordinate
    - id: unmap_contaminant_reads
      label: Unmap contaminant reads
      doc: |-
        Detect reads originating from foreign organisms (e.g. bacterial DNA in a non-bacterial sample), and unmap + label those reads accordingly.
      type: boolean?
      inputBinding:
        prefix: --UNMAP_CONTAMINANT_READS
        position: 4
        shellQuote: false
      sbg:altPrefix: -UNMAP_CONTAM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: unmapped_bam
      label: Unmapped BAM/SAM file
      doc: Original SAM or BAM file of unmapped reads, which must be in queryname
        order.
      type: File
      inputBinding:
        prefix: --UNMAPPED_BAM
        position: 4
        shellQuote: false
      sbg:altPrefix: -UNMAPPED
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM
    - id: unmapped_read_strategy
      label: Unmapped read strategy
      doc: |-
        How to deal with alignment information in reads that are being unmapped (e.g. due to cross-species contamination.) Currently ignored unless UNMAP_CONTAMINANT_READS = true
      type:
      - 'null'
      - name: unmapped_read_strategy
        type: enum
        symbols:
        - COPY_TO_TAG
        - DO_NOT_CHANGE
        - MOVE_TO_TAG
      inputBinding:
        prefix: --UNMAPPED_READ_STRATEGY
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: DO_NOT_CHANGE
    - id: validation_stringency
      label: Validation stringency
      doc: |-
        Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
      type:
      - 'null'
      - name: validation_stringency
        type: enum
        symbols:
        - STRICT
        - LENIENT
        - SILENT
      inputBinding:
        prefix: --VALIDATION_STRINGENCY
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
    - id: output_prefix
      label: Output prefix
      doc: Output file name prefix.
      type: string?
      sbg:category: Optional Parameters
    - id: output_file_format
      label: Output file format
      doc: Output file format
      type:
      - 'null'
      - name: output_file_format
        type: enum
        symbols:
        - bam
        - sam
      sbg:category: Optional parameters
    - id: cpu_per_job
      label: CPU per job
      doc: CPU per job.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'

    outputs:
    - id: out_alignments
      label: Output merged SAM or BAM file
      doc: Output merged SAM or BAM file.
      type: File
      secondaryFiles:
      - |-
        ${
            if (self.nameext == ".bam" && inputs.create_index)
            {
                return [self.basename + ".bai", self.nameroot + ".bai"];
            }
            else {
                return []; 
            }
        }
      outputBinding:
        glob: '*am'
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: SAM, BAM

    baseCommand: []
    arguments:
    - prefix: ''
      position: 0
      valueFrom: /opt/gatk
      shellQuote: false
    - prefix: ''
      position: 1
      valueFrom: |-
        ${
            if (inputs.memory_per_job)
            {
                return "--java-options";
            }
            else {
                return '';
            }
        }
      shellQuote: false
    - prefix: ''
      position: 2
      valueFrom: |-
        ${
            if (inputs.memory_per_job) {
                return '\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\"';
            }
            else {
                return ''; 
                
            }
        }
      shellQuote: false
    - position: 3
      valueFrom: MergeBamAlignment
      shellQuote: false
    - prefix: ''
      position: 4
      valueFrom: |-
        ${
            var in_alignments = [].concat(inputs.in_alignments);
            var output_ext = inputs.output_file_format ? inputs.output_file_format : in_alignments[0].path.split('.').pop();
            var output_prefix = '';
            var file1_name = ''; 
            var file2_name = ''; 
            if (inputs.output_prefix)
            {
                output_prefix = inputs.output_prefix;
            }
            else 
            {
                if (in_alignments.length > 1)
                {
                    in_alignments.sort(function(file1, file2) {
                        file1_name = file1.path.split('/').pop().toUpperCase();
                        file2_name = file2.path.split('/').pop().toUpperCase();
                        if (file1_name < file2_name) {
                            return -1;
                        }
                        if (file1_name > file2_name) {
                            return 1;
                        }
                        // names must be equal
                        return 0;
                    });
                }
                
                var in_alignments_first =  in_alignments[0];
                if (in_alignments_first.metadata && in_alignments_first.metadata.sample_id)
                {
                    output_prefix = in_alignments_first.metadata.sample_id;
                }
                else 
                {
                    output_prefix = in_alignments_first.path.split('/').pop().split('.')[0];
                }
                
                if (in_alignments.length > 1)
                {
                    output_prefix = output_prefix + "." + in_alignments.length;
                }
            }
            
            return "--OUTPUT " + output_prefix + ".merged." + output_ext;
        }
      shellQuote: false
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-mergebamalignment-4-1-0-0/14
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    sbg:content_hash: a758b43167e957642f45a0aad07716ff3b8c8c6a379cf76b35f10b0a3f5a121b8
    sbg:contributors:
    - uros_sipetic
    - nemanja.vucic
    - nens
    - veliborka_josipovic
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/37
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552666475
    sbg:id: h-0ceca83e/h-b97d0632/h-8d60708c/0
    sbg:image_url:
    sbg:latestRevision: 14
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - id: https://software.broadinstitute.org/gatk/
      label: Homepage
    - id: https://github.com/broadinstitute/gatk/
      label: Source Code
    - id: |-
        https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
      label: Download
    - id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
      label: Publications
    - id: |-
        https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_MergeSamFiles.php
      label: Documentation
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1560336165
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 14
    sbg:revisionNotes: |-
      Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/37
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1552666475
      sbg:revision: 0
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/12
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554492767
      sbg:revision: 1
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/23
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554720890
      sbg:revision: 2
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/24
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554999266
      sbg:revision: 3
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/25
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557734540
      sbg:revision: 4
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/26
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000585
      sbg:revision: 5
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/27
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558017849
      sbg:revision: 6
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/28
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351570
      sbg:revision: 7
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/29
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558370509
      sbg:revision: 8
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/30
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558427482
      sbg:revision: 9
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/31
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558448356
      sbg:revision: 10
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/32
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558453788
      sbg:revision: 11
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/33
    - sbg:modifiedBy: nemanja.vucic
      sbg:modifiedOn: 1559750464
      sbg:revision: 12
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/34
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1560335266
      sbg:revision: 13
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/36
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1560336165
      sbg:revision: 14
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/37
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: out_alignments
  sbg:x: -9
  sbg:y: 53.96965026855469
- id: gatk_samtofastq_4_1_0_0
  label: GATK SamToFastq
  in:
  - id: include_non_pf_reads
    default: true
  - id: in_alignments
    source: in_alignments
  - id: interleave
    default: true
  scatter:
  - in_alignments
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK SamToFastq
    doc: |-
      The **GATK SamToFastq** tool converts a SAM or BAM file to FASTQ.

      This tool extracts read sequences and qualities from the input SAM/BAM file and writes them into the output file in Sanger FASTQ format.

      In the RC mode (default is True), if the read is aligned and the alignment is to the reverse strand on the genome, the read sequence from input SAM file will be reverse-complemented prior to writing it to FASTQ in order to correctly restore the original read sequence as it was generated by the sequencer [1].

      *A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*

      ###Common Use Cases

      * The **GATK SamToFastq** tool requires a BAM/SAM file on its **Input BAM/SAM file** (`--INPUT`) input. The tool generates a single-end FASTQ file on its **Output FASTQ file(s)** output if the input BAM/SAM file is single end. In case the input file is paired end, the tool outputs the first end of the pair FASTQ and the second end of the pair FASTQ on its **Output FASTQ file(s)** output, except when the **Interleave** (`--INTERLEAVE`) option is set to True. If the output is an interleaved FASTQ file, if paired, each line will have /1 or /2 to describe which end it came from.

      * The **GATK SamToFastq** tool supports an optional parameter  **Output by readgroup** (`--OUTPUT_BY_READGROUP`) which, when true, outputs a FASTQ file per read group (two FASTQ files per read group if the group is paired).

      * Usage example (input BAM file is single-end):

      ```
      gatk SamToFastq 
           --INPUT input.bam
           --FASTQ output.fastq
      ```





      * Usage example (input BAM file is paired-end):

      ```
      gatk SamToFastq 
           --INPUT input.bam
           --FASTQ output.pe_1.fastq
           --SECOND_END_FASTQ output.pe_2.fastq
           --UNPAIRED_FASTQ unpaired.fastq

      ```

      ###Changes Introduced by Seven Bridges

      * The GATK SamToFastq tool is implemented to check if the input alignments file contains single-end or paired-end data and according to that generates different command lines for these two modes and thus produces appropriate output files on its **Output FASTQ file(s)** output (one FASTQ file in single-end mode and two FASTQ files if the input alignment file contains paired-end data). 

      * All output files will be prefixed using the **Output prefix** parameter. In case the **Output prefix** is not provided, the output prefix will be the same as the Sample ID metadata from the **input SAM/BAM file**, if the Sample ID metadata exists. Otherwise, the output prefix will be inferred from the **Input SAM/BAM** filename. This way, having identical names of the output files between runs is avoided.

      * For paired-end read files, in order to successfully run alignment with STAR, this tool adds the appropriate **paired-end** metadata field in the output FASTQ files.

      ###Common Issues and Important Notes

      * None

      ###Performance Benchmarking

      Below is a table describing runtimes and task costs of **GATK SamToFastq** for a couple of different samples, executed on the AWS cloud instances:

      | Experiment type |  Input size | Paired-end | # of reads | Read length | Duration |  Cost | Instance (AWS) | 
      |:--------------:|:------------:|:--------:|:-------:|:---------:|:----------:|:------:|:------:|
      |     RNA-Seq     |  1.9 GB |     Yes    |     16M     |     101     |   4min   | ~0.03$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     |  5.7 GB |     Yes    |     50M     |     101     |   7min   | ~0.04$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 9.5 GB |     Yes    |     82M    |     101     |  10min  | ~0.07$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 19 GB |     Yes    |     164M    |     101     |  20min  | ~0.13$ | c4.2xlarge (8 CPUs) |

      *Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*


      ###References

      [1] [GATK SamToFastq](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_sam_SamToFastq)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: "${\n    return inputs.cpu_per_job ? inputs.cpu_per_job : 1;\n}"
      ramMin: |-
        ${
            var memory = 4096;
            if (inputs.memory_per_job) 
            {
                memory = inputs.memory_per_job;
            }
            if (inputs.memory_overhead_per_job)
            {
                memory += inputs.memory_overhead_per_job;
            }
            return memory;
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file)) {
                file['metadata'] = {}
            }
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

    inputs:
    - id: clipping_action
      label: Clipping action
      doc: |-
        The action that should be taken with clipped reads: 'X' means the reads and qualities should be trimmed at the clipped position; 'N' means the bases should be changed to Ns in the clipped region; and any integer means that the base qualities should be set to that value in the clipped region.
      type: string?
      inputBinding:
        prefix: --CLIPPING_ACTION
        position: 5
        shellQuote: false
      sbg:altPrefix: -CLIP_ACT
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: clipping_attribute
      label: Clipping attribute
      doc: |-
        The attribute that stores the position at which the SAM record should be clipped.
      type: string?
      inputBinding:
        prefix: --CLIPPING_ATTRIBUTE
        position: 5
        shellQuote: false
      sbg:altPrefix: -CLIP_ATTR
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: clipping_min_length
      label: Clipping min length
      doc: |-
        When performing clipping with the CLIPPING_ATTRIBUTE and CLIPPING_ACTION parameters, ensure that the resulting reads after clipping are at least CLIPPING_MIN_LENGTH bases long. If the original read is shorter than CLIPPING_MIN_LENGTH then the original read length will be maintained.
      type: int?
      inputBinding:
        prefix: --CLIPPING_MIN_LENGTH
        position: 5
        shellQuote: false
      sbg:altPrefix: -CLIP_MIN
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: compress_outputs_per_rg
      label: Compress outputs per RG
      doc: |-
        Compress output FASTQ files per read group using gzip and append a .gz extension to the file names. Cannot be used in conjuction with argument(s) FASTQ (F) SECOND_END_FASTQ (F2) UNPAIRED_FASTQ (FU).
      type: boolean?
      inputBinding:
        prefix: --COMPRESS_OUTPUTS_PER_RG
        position: 5
        shellQuote: false
      sbg:altPrefix: -GZOPRG
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: compression_level
      label: Compression level
      doc: Compression level for all compressed files created (e.g. BAM and VCF).
      type: int?
      inputBinding:
        prefix: --COMPRESSION_LEVEL
        position: 5
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
    - id: include_non_pf_reads
      label: Include non PF reads
      doc: |-
        Include non-PF reads from the SAM file into the output FASTQ files. PF means 'passes filtering'. Reads whose 'not passing quality controls' flag is set are non-PF reads. See GATK Dictionary for more info.
      type: boolean?
      inputBinding:
        prefix: --INCLUDE_NON_PF_READS
        position: 5
        shellQuote: false
      sbg:altPrefix: -NON_PF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: include_non_primary_alignments
      label: Include non primary alignments
      doc: |-
        If true, include non-primary alignments in the output. Support of non-primary alignments in SamToFastq is not comprehensive, so there may be exceptions if this is set to true and there are paired reads with non-primary alignments.
      type: boolean?
      inputBinding:
        prefix: --INCLUDE_NON_PRIMARY_ALIGNMENTS
        position: 5
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: in_alignments
      label: Input SAM/BAM file
      doc: Input SAM/BAM file to extract reads from.
      type: File
      inputBinding:
        prefix: --INPUT
        position: 5
        shellQuote: false
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: SAM, BAM
    - id: interleave
      label: Interleave
      doc: |-
        Will generate an interleaved FASTQ if paired, each line will have /1 or /2 to describe which end it came from.
      type: boolean?
      inputBinding:
        prefix: --INTERLEAVE
        position: 5
        shellQuote: false
      sbg:altPrefix: -INTER
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: memory_overhead_per_job
      label: Memory overhead per job
      doc: |-
        This input allows a user to set the desired overhead memory when running a tool or adding it to a workflow. This amount will be added to the Memory per job in the Memory requirements section but it will not be added to the -Xmx parameter leaving some memory not occupied which can be used as stack memory (-Xmx parameter defines heap memory). This input should be defined in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      type: int?
      sbg:category: Platform Options
    - id: memory_per_job
      label: Memory per job
      doc: |-
        This input allows a user to set the desired memory requirement when running a tool or adding it to a workflow. This value should be propagated to the -Xmx parameter too.This input should be defined in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: 2048 MB
    - id: output_per_rg
      label: Output per RG
      doc: |-
        Output a FASTQ file per read group (two FASTQ files per read group if the group is paired). Cannot be used in conjuction with argument(s)FASTQ (F) SECOND_END_FASTQ (F2) UNPAIRED_FASTQ (FU).
      type: boolean?
      inputBinding:
        prefix: --OUTPUT_PER_RG
        position: 5
        shellQuote: false
      sbg:altPrefix: -OPRG
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: quality
      label: Quality
      doc: End-trim reads using the phred/bwa quality trimming algorithm and this
        quality.
      type: int?
      inputBinding:
        prefix: --QUALITY
        position: 5
        shellQuote: false
      sbg:altPrefix: -Q
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: re_reverse
      label: Re reverse
      doc: |-
        Re-reverse bases and qualities of reads with negative strand flag set before writing them to FASTQ.
      type:
      - 'null'
      - name: re_reverse
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --RE_REVERSE
        position: 5
        shellQuote: false
      sbg:altPrefix: -RC
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: read1_max_bases_to_write
      label: Read1 max bases to write
      doc: |-
        The maximum number of bases to write from read 1 after trimming. If there are fewer than this many bases left after trimming, all will be written. If this value is null then all bases left after trimming will be written.
      type: int?
      inputBinding:
        prefix: --READ1_MAX_BASES_TO_WRITE
        position: 5
        shellQuote: false
      sbg:altPrefix: -R1_MAX_BASES
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: read1_trim
      label: Read1 trim
      doc: The number of bases to trim from the beginning of read 1.
      type: int?
      inputBinding:
        prefix: --READ1_TRIM
        position: 5
        shellQuote: false
      sbg:altPrefix: -R1_TRIM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: read2_max_bases_to_write
      label: Read2 max bases to write
      doc: |-
        The maximum number of bases to write from read 2 after trimming. If there are fewer than this many bases left after trimming, all will be written. If this value is null then all bases left after trimming will be written.
      type: int?
      inputBinding:
        prefix: --READ2_MAX_BASES_TO_WRITE
        position: 5
        shellQuote: false
      sbg:altPrefix: -R2_MAX_BASES
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: read2_trim
      label: Read2 trim
      doc: The number of bases to trim from the beginning of read 2.
      type: int?
      inputBinding:
        prefix: --READ2_TRIM
        position: 5
        shellQuote: false
      sbg:altPrefix: -R2_TRIM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: rg_tag
      label: RG tag
      doc: The read group tag (PU or ID) to be used to output a FASTQ file per read
        group.
      type: string?
      inputBinding:
        prefix: --RG_TAG
        position: 5
        shellQuote: false
      sbg:altPrefix: -RGT
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: PU
    - id: validation_stringency
      label: Validation stringency
      doc: |-
        Validation stringency for all SAM files read by this program. Setting stringency to silent can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
      type:
      - 'null'
      - name: validation_stringency
        type: enum
        symbols:
        - STRICT
        - LENIENT
        - SILENT
      inputBinding:
        prefix: --VALIDATION_STRINGENCY
        position: 5
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
    - id: output_prefix
      label: Output prefix
      doc: Output file name prefix.
      type: string?
      sbg:category: Optional Arguments
    - id: compress_outputs
      label: Compress output file(s)
      doc: Compress output file(s).
      type: boolean?
      sbg:category: Optional parameters
      sbg:toolDefaultValue: 'false'
    - id: cpu_per_job
      label: CPU per job
      doc: CPU per job.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'

    outputs:
    - id: out_reads
      label: Output FASTQ file(s)
      doc: Output FASTQ file(s).
      type: File[]?
      outputBinding:
        glob: |-
          ${
              var output_ext = inputs.compress_outputs ? ".fastq.gz" : ".fastq";
              var interleave = inputs.interleave;
              if (!inputs.outputs_by_readgroup)
              {
                  if (interleave)
                  {
                      return "*interleaved" + output_ext;
                  }
                  else
                  {
                      return ["*pe_1" + output_ext, "*pe_2" + output_ext, "*se" + output_ext];
                  }

              }
              else
              {
                  return "*" + output_ext;
              }
          }
        outputEval: |-
          ${ 
              self = [].concat(self)
              
              function getPairedEnd(filename)
              {
                  if (filename.lastIndexOf(".fastq") !== 0 && filename[filename.lastIndexOf(".fastq") - 2 ]=="_") 
                  {
                      return filename[filename.lastIndexOf(".fastq") - 1 ];
                  } 
                  else 
                  {
                      return "";
                  }
              }
              
              var out = inheritMetadata(self,inputs.in_alignments);
              for (var i=0; i < out.length; i++)
              {
                  out[i].metadata['paired_end'] = getPairedEnd(out[i].path);
              }
              
              return out;
          }
      sbg:fileTypes: FASTQ, FASTQ.GZ
    - id: unmapped_reads
      label: Unpaired reads
      doc: Unpaired reads.
      type: File?
      outputBinding:
        glob: |-
          ${
              var output_ext = inputs.compress_outputs ? ".fastq.gz" : ".fastq";
              var interleave = inputs.interleave;
              if (!inputs.outputs_by_readgroup)
              {
                  if (!interleave)
                  {
                      return "*unpaired" + output_ext;
                  }
                  else 
                  {
                       return ""; 
                  }      
              }
            else {
                 return "";  
            
            }
          }
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: FASTQ, FASTQ.GZ

    baseCommand: []
    arguments:
    - prefix: ''
      position: 0
      valueFrom: |+
        ${
            var in_alignments = [].concat(inputs.in_alignments)[0];
            var output_ext    = inputs.compress_outputs ? ".fastq.gz" : ".fastq";
            var interleave    = inputs.interleave;
            var output_prefix = ''; 
            var cmd_line      = '';
            cmd_line          = "cmd='' && paired_end=`samtools view -h " + in_alignments.path + " | head -n 500000 | samtools view -Sc -f 0x1 -`";

          if (!inputs.outputs_by_readgroup)
            {
                if (inputs.output_prefix)
                {
                    output_prefix = inputs.output_prefix;
                }
                else
                {
                    if (in_alignments.metadata && in_alignments.metadata.sample_id)
                    {
                        output_prefix = in_alignments.metadata.sample_id;
                    }
                    else
                    {
                        output_prefix = in_alignments.path.split('/').pop().split('.')[0];
                    }          
                }           
                
                cmd_line = cmd_line + " && if [ $paired_end != 0 ]; then cmd='--FASTQ " + output_prefix; 
                
                if (interleave)
                {
                    cmd_line = cmd_line + ".interleaved" + output_ext + "';";
                }
                else
                {
                    cmd_line = cmd_line + ".pe_1" + output_ext;
                    cmd_line = cmd_line + " --SECOND_END_FASTQ " + output_prefix + ".pe_2" + output_ext;
                    cmd_line = cmd_line + " --UNPAIRED_FASTQ " + output_prefix + ".unpaired" + output_ext + "';";
                }        
                cmd_line = cmd_line + " else cmd='--FASTQ " + output_prefix  + ".se" + output_ext + "'; fi;";
                return cmd_line;
            }
            else
            {
                return "cmd='--OUTPUT_DIR .'";
            }
        }

      shellQuote: false
    - prefix: ''
      position: 1
      valueFrom: /opt/gatk
      shellQuote: false
    - prefix: ''
      position: 2
      valueFrom: |-
        ${
            if (inputs.memory_per_job)
            {
                return "--java-options";
            }
            else {
                return '';
            }
            
        }
      shellQuote: false
    - prefix: ''
      position: 3
      valueFrom: |-
        ${
            if (inputs.memory_per_job) {
                return '\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\"';
            }
            else {
                return "";
            }
        }
      shellQuote: false
    - prefix: ''
      position: 4
      valueFrom: SamToFastq
      shellQuote: false
    - prefix: ''
      position: 6
      valueFrom: "${\n        return '$cmd';\n}"
      shellQuote: false
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-samtofastq-4-1-0-0/15
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    sbg:content_hash: a21e1194e724a1f17bceabd4d2040324713c2a5c63896adcebbc777578b2bfef5
    sbg:contributors:
    - uros_sipetic
    - nens
    - veliborka_josipovic
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552663400
    sbg:id: h-5c321f4c/h-45f450f4/h-436e400a/0
    sbg:image_url:
    sbg:latestRevision: 15
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - id: https://software.broadinstitute.org/gatk/
      label: Homepage
    - id: https://github.com/broadinstitute/gatk/
      label: Source Code
    - id: |-
        https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
      label: Download
    - id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
      label: Publications
    - id: |-
        https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_SamToFastq.php
      label: Documentation
    sbg:modifiedBy: veliborka_josipovic
    sbg:modifiedOn: 1561548030
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 15
    sbg:revisionNotes: Added glob for single end output fastq
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1552663400
      sbg:revision: 0
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/19
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1552663734
      sbg:revision: 1
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/20
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554492676
      sbg:revision: 2
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/27
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554493243
      sbg:revision: 3
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/28
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554720826
      sbg:revision: 4
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/29
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554999298
      sbg:revision: 5
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/30
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557484228
      sbg:revision: 6
      sbg:revisionNotes: Updated Description
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557745933
      sbg:revision: 7
      sbg:revisionNotes: ''
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557918579
      sbg:revision: 8
      sbg:revisionNotes: 'v32: [input]'
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557919927
      sbg:revision: 9
      sbg:revisionNotes: v5->update
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558002424
      sbg:revision: 10
      sbg:revisionNotes: output required
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558015833
      sbg:revision: 11
      sbg:revisionNotes: unmapped_reads required
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558022954
      sbg:revision: 12
      sbg:revisionNotes: strict js for glob
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558023482
      sbg:revision: 13
      sbg:revisionNotes: strict javascript for unmapped_reads
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558354170
      sbg:revision: 14
      sbg:revisionNotes: return '';
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1561548030
      sbg:revision: 15
      sbg:revisionNotes: Added glob for single end output fastq
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: out_reads
  - id: unmapped_reads
  sbg:x: -444.0947265625
  sbg:y: 120.06857299804688
- id: gatk_sortsam_4_1_0_0
  label: GATK SortSam
  in:
  - id: in_alignments
    source: gatk_markduplicates_4_1_0_0/out_alignments
  - id: sort_order
    default: coordinate
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK SortSam
    doc: |-
      The **GATK SortSam** tool sorts the input SAM or BAM file by coordinate, queryname (QNAME), or some other property of the SAM record.

      The **GATK SortOrder** of a SAM/BAM file is found in the SAM file header tag @HD in the field labeled SO.  For a coordinate
      sorted SAM/BAM file, read alignments are sorted first by the reference sequence name (RNAME) field using the reference
      sequence dictionary (@SQ tag).  Alignments within these subgroups are secondarily sorted using the left-most mapping
      position of the read (POS).  Subsequent to this sorting scheme, alignments are listed arbitrarily.</p><p>For
      queryname-sorted alignments, all alignments are grouped using the queryname field but the alignments are not necessarily
      sorted within these groups.  Reads having the same queryname are derived from the same template


      ###Common Use Cases

      The **GATK SortSam** tool requires a BAM/SAM file on its **Input SAM/BAM file**   (`--INPUT`)  input. The tool sorts input file in the order defined by (`--SORT_ORDER`) parameter. Available sort order options are `queryname`, `coordinate` and `duplicate`.  

      * Usage example:

      ```
      java -jar picard.jar SortSam
           --INPUT=input.bam 
           --SORT_ORDER=coordinate
      ```


      ###Changes Introduced by Seven Bridges

      * Prefix of the output file is defined with the optional parameter **Output prefix**. If **Output prefix** is not provided, name of the sorted file is obtained from **Sample ID** metadata from the **Input SAM/BAM file**, if the **Sample ID** metadata exists. Otherwise, the output prefix will be inferred form the **Input SAM/BAM file** filename. 


      ###Common Issues and Important Notes

      * None


      ###Performance Benchmarking
      Below is a table describing runtimes and task costs of **GATK SortSam** for a couple of different samples, executed on the AWS cloud instances:

      | Experiment type |  Input size | Paired-end | # of reads | Read length | Duration |  Cost | Instance (AWS) | 
      |:--------------:|:------------:|:--------:|:-------:|:---------:|:----------:|:------:|:------:|
      |     WGS     |          |     Yes    |     16M     |     101     |   4min   | ~0.03$ | c4.2xlarge (8 CPUs) | 
      |     WGS     |         |     Yes    |     50M     |     101     |   7min   | ~0.04$ | c4.2xlarge (8 CPUs) | 
      |     WGS     |         |     Yes    |     82M    |     101     |  10min  | ~0.07$ | c4.2xlarge (8 CPUs) | 
      |     WES     |         |     Yes    |     164M    |     101     |  20min  | ~0.13$ | c4.2xlarge (8 CPUs) |

      *Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*



      ###References
      [1] [GATK SortSam home page](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_sam_SortSam.php)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: "${\n    return inputs.cpu_per_job ? inputs.cpu_per_job : 1;\n}"
      ramMin: |-
        ${
            var memory = 4096;
            if (inputs.memory_per_job) 
            {
                memory = inputs.memory_per_job;
            }
            if (inputs.memory_overhead_per_job)
            {
                memory += inputs.memory_overhead_per_job;
            }
            return memory;
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


        var setMetadata = function(file, metadata) {
            if (!('metadata' in file)) {
                file['metadata'] = {}
            }
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

    inputs:
    - id: in_alignments
      label: Input SAM/BAM file
      doc: Input BAM or SAM file to sort.  Required
      type: File
      inputBinding:
        prefix: --INPUT
        position: 4
        shellQuote: false
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM
    - id: output_prefix
      label: Output prefix
      doc: Sorted bam or sam output file.
      type: string?
      sbg:altPrefix: -O
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: sample_id.sorted.bam
    - id: compression_level
      label: Compression level
      doc: Compression level for all compressed files created (e.g. Bam and vcf).
      type: int?
      inputBinding:
        prefix: --COMPRESSION_LEVEL
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
    - id: create_index
      label: Create index
      doc: Whether to create a bam index when writing a coordinate-sorted bam file.
      type: boolean?
      inputBinding:
        prefix: --CREATE_INDEX
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: create_md5_file
      label: Create md5 file
      doc: Whether to create an md5 digest for any bam or fastq files created.
      type: boolean?
      inputBinding:
        prefix: --CREATE_MD5_FILE
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: max_records_in_ram
      label: Max records in ram
      doc: |-
        When writing files that need to be sorted, this will specify the number of records stored in ram before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of ram needed.
      type: int?
      inputBinding:
        prefix: --MAX_RECORDS_IN_RAM
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
    - id: validation_stringency
      label: Validation stringency
      doc: |-
        Validation stringency for all sam files read by this program. Setting stringency to silent can improve performance when processing a bam file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
      type:
      - 'null'
      - name: validation_stringency
        type: enum
        symbols:
        - STRICT
        - LENIENT
        - SILENT
      inputBinding:
        prefix: --VALIDATION_STRINGENCY
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
    - id: memory_per_job
      label: Memory Per Job
      doc: Memory which will be allocated for execution.
      type: int?
      sbg:category: Execution
    - id: memory_overhead_per_job
      label: Memory Overhead Per Job
      doc: Memory overhead which will be allocated for one job.
      type: int?
      sbg:category: Execution
    - id: sort_order
      doc: |-
        Sort order of output file.   Required. Possible values: {
                                      queryname (Sorts according to the readname. This will place read-pairs and other derived
                                      reads (secondary and supplementary) adjacent to each other. Note that the readnames are
                                      compared lexicographically, even though they may include numbers. In paired reads, Read1
                                      sorts before Read2.)
                                      coordinate (Sorts primarily according to the SEQ and POS fields of the record. The
                                      sequence will sorted according to the order in the sequence dictionary, taken from from
                                      the header of the file. Within each reference sequence, the reads are sorted by the
                                      position. Unmapped reads whose mates are mapped will be placed near their mates. Unmapped
                                      read-pairs are placed after all the mapped reads and their mates.)
                                      duplicate (Sorts the reads so that duplicates reads are adjacent. Required that the
                                      mate-cigar (MC) tag is present. The resulting will be sorted by library, unclipped 5-prime
                                      position, orientation, and mate's unclipped 5-prime position.)
                                      }
      type:
        name: sort_order
        type: enum
        symbols:
        - queryname
        - coordinate
        - duplicate
      inputBinding:
        prefix: --SORT_ORDER
        position: 7
        shellQuote: false
      sbg:altPrefix: -SO
      sbg:category: Required  Arguments
    - id: cpu_per_job
      label: CPU per job
      doc: |-
        This input allows a user to set the desired CPU requirement when running a tool or adding it to a workflow.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
    - id: output_file_format
      label: Output file format
      doc: Output file format.
      type:
      - 'null'
      - name: output_file_format
        type: enum
        symbols:
        - bam
        - sam
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: Same as input

    outputs:
    - id: out_alignments
      label: Sorted BAM/SAM
      doc: Sorted BAM or SAM output file.
      type: File?
      secondaryFiles:
      - |-
        ${
           if (inputs.create_index)
           {
               return [self.basename + ".bai", self.nameroot + ".bai"]
           }
           else {
               return []; 
           }
        }
      outputBinding:
        glob: '*am'
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: BAM, SAM

    baseCommand: []
    arguments:
    - position: 0
      valueFrom: /opt/gatk
      shellQuote: false
    - position: 1
      valueFrom: --java-options
      shellQuote: false
    - prefix: ''
      position: 2
      valueFrom: |-
        ${
            if (inputs.memory_per_job) {
                return '\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\"';
            }
            return '\"-Xmx2048M\"';
        }
      shellQuote: false
    - position: 3
      valueFrom: SortSam
      shellQuote: false
    - prefix: ''
      position: 4
      valueFrom: |-
        ${
            var tmp = [].concat(inputs.in_alignments);
            var ext = '';
          
            if (inputs.output_file_format){
                ext = inputs.output_file_format;
            }    else {
                ext = tmp[0].path.split(".").pop();
            }
            
            
            if (inputs.output_prefix) {
                return '-O ' +  inputs.output_prefix + ".sorted." + ext;
              
            }else if (tmp[0].metadata && tmp[0].metadata.sample_id) {
                
                return '-O ' +  tmp[0].metadata.sample_id + ".sorted." + ext;
            } else {
                 
                return '-O ' +  tmp[0].path.split('/').pop().split(".")[0] + ".sorted."+ext;
            }
            
            
        }
      shellQuote: false
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-sortsam-4-1-0-0/8
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    sbg:content_hash: a4d21247730823bddd1b0c24a25cc7b27bea6e061eacc901c23e642f333f458d5
    sbg:contributors:
    - nens
    - uros_sipetic
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/19
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1555498331
    sbg:id: h-4c2f17d9/h-9140f101/h-f3a9ebb3/0
    sbg:image_url:
    sbg:latestRevision: 8
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - id: https://software.broadinstitute.org/gatk/
      label: Homepage
    - id: |-
        https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_SortSam.php
      label: Documentation
    - id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
      label: Publications
    - id: https://github.com/broadinstitute/gatk/
      label: Source code
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1561632457
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 8
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/19
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1555498331
      sbg:revision: 0
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/2
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1555582270
      sbg:revision: 1
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/9
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557417459
      sbg:revision: 2
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/11
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557734528
      sbg:revision: 3
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/13
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000570
      sbg:revision: 4
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/14
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558009951
      sbg:revision: 5
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/15
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351565
      sbg:revision: 6
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/17
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558449641
      sbg:revision: 7
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/18
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1561632457
      sbg:revision: 8
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/19
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: out_alignments
  sbg:x: 434.41656494140625
  sbg:y: 186.55223083496094
- id: gatk_setnmmdanduqtags_4_1_0_0
  label: GATK SetNmMdAndUqTags
  in:
  - id: create_index
    default: true
  - id: in_alignments
    source: gatk_sortsam_4_1_0_0/out_alignments
  - id: reference_sequence
    source: in_reference
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK SetNmMdAndUqTags
    doc: |-
      The **GATK SetNmMdAndUqTags** tool takes in a coordinate-sorted SAM or BAM and calculatesthe NM, MD, and UQ tags by comparing it with the reference. 

      The **GATK SetNmMdAndUqTags**  may be needed when **GATK MergeBamAlignment** was run with **SORT_ORDER** other than `coordinate` and thus could not fix these tags. 


      ###Common Use Cases
      The **GATK SetNmMdAndUqTags** tool  fixes NM, MD and UQ tags in SAM/BAM file **Input SAM/BAM file**   (`--INPUT`)  input. This tool takes in a coordinate-sorted SAM or BAM file and calculates the NM, MD, and UQ tags by comparing with the reference **Reference sequence** (`--REFERENCE_SEQUENCE`).

      * Usage example:

      ```
      java -jar picard.jar SetNmMdAndUqTags
           --REFERENCE_SEQUENCE=reference_sequence.fasta
           --INPUT=sorted.bam
      ```


      ###Changes Introduced by Seven Bridges

      * Prefix of the output file is defined with the optional parameter **Output prefix**. If **Output prefix** is not provided, name of the sorted file is obtained from **Sample ID** metadata form the **Input SAM/BAM file**, if the **Sample ID** metadata exists. Otherwise, the output prefix will be inferred form the **Input SAM/BAM file** filename. 



      ###Common Issues and Important Notes

      * The **Input SAM/BAM file** must be coordinate sorted in order to run  **GATK SetNmMdAndUqTags**. 
      * If specified, the MD and NM tags can be ignored and only the UQ tag be set. 


      ###References
      [1] [GATK SetNmMdAndUqTags home page](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_sam_SetNmMdAndUqTags.php)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: "${\n    return inputs.cpu_per_job ? inputs.cpu_per_job : 1;\n}"
      ramMin: |-
        ${
            var memory = 4096;
            if (inputs.memory_per_job) 
            {
                memory = inputs.memory_per_job;
            }
            if (inputs.memory_overhead_per_job)
            {
                memory += inputs.memory_overhead_per_job;
            }
            return memory;
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


        var setMetadata = function(file, metadata) {
            if (!('metadata' in file)) {
                file['metadata'] = {}
            }
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

    inputs:
    - id: validation_stringency
      label: Validation stringency
      doc: |-
        Validation stringency for all sam files read by this program. Setting stringency to silent can improve performance when processing a bam file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
      type:
      - 'null'
      - name: validation_stringency
        type: enum
        symbols:
        - STRICT
        - LENIENT
        - SILENT
      inputBinding:
        prefix: --VALIDATION_STRINGENCY
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
    - id: output_prefix
      label: Output
      doc: The fixed bam or sam output prefix name.
      type: string?
      sbg:altPrefix: -O
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: sample_id.fixed.bam
    - id: memory_overhead_per_job
      label: Memory Overhead Per Job
      doc: |-
        This input allows a user to set the desired overhead memory when running a tool or adding it to a workflow. This amount will be added to the Memory per job in the Memory requirements section but it will not be added to the -Xmx parameter leaving some memory not occupied which can be used as stack memory (-Xmx parameter defines heap memory). This input should be defined in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      type: int?
      sbg:category: Execution
    - id: max_records_in_ram
      label: Max records in ram
      doc: |-
        When writing files that need to be sorted, this will specify the number of records stored in ram before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of ram needed.
      type: int?
      inputBinding:
        prefix: --MAX_RECORDS_IN_RAM
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
    - id: create_index
      label: Create index
      doc: Whether to create a bam index when writing a coordinate-sorted bam file.
      type: boolean?
      inputBinding:
        prefix: --CREATE_INDEX
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: is_bisulfite_sequence
      label: Is bisulfite sequence
      doc: Whether the file contains bisulfite sequence (used when calculating the
        nm tag).
      type: boolean?
      inputBinding:
        prefix: --IS_BISULFITE_SEQUENCE
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: compression_level
      label: Compression level
      doc: Compression level for all compressed files created (e.g. Bam and vcf).
      type: int?
      inputBinding:
        prefix: --COMPRESSION_LEVEL
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
    - id: memory_per_job
      label: Memory Per Job
      doc: |-
        This input allows a user to set the desired memory requirement when running a tool or adding it to a workflow. This value should be propagated to the -Xmx parameter too.This input should be defined in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      type: int?
      sbg:category: Execution
    - id: in_alignments
      label: Input SAM/BAM file
      doc: The BAM or SAM file to fix.
      type: File
      inputBinding:
        prefix: --INPUT
        position: 4
        shellQuote: false
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM
    - id: reference_sequence
      label: Reference sequence
      doc: Reference sequence FASTA file.
      type: File
      inputBinding:
        prefix: --REFERENCE_SEQUENCE
        position: 4
        shellQuote: false
      sbg:altPrefix: -R
      sbg:category: Required Arguments
      sbg:fileTypes: FASTA, FA
    - id: set_only_uq
      label: Set only uq
      doc: Only set the uq tag, ignore md and nm.
      type: boolean?
      inputBinding:
        prefix: --SET_ONLY_UQ
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: cpu_per_job
      label: CPU per job
      doc: |-
        This input allows a user to set the desired CPU requirement when running a tool or adding it to a workflow.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
    - id: output_file_format
      label: Output file format
      doc: Output file format.
      type:
      - 'null'
      - name: output_file_format
        type: enum
        symbols:
        - bam
        - sam
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: Same as input

    outputs:
    - id: out_alignments
      label: Output BAM/SAM file
      doc: Output BAM/SAM file with fixed tags.
      type: File?
      secondaryFiles:
      - |-
        ${  
            if (inputs.create_index)
            {
                return self.nameroot + ".bai";
            }
            else {
                return ''; 
            }
        }
      outputBinding:
        glob: '*am'
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: BAM, SAM

    baseCommand: []
    arguments:
    - position: 0
      valueFrom: /opt/gatk
      shellQuote: false
    - position: 1
      valueFrom: --java-options
      shellQuote: false
    - prefix: ''
      position: 2
      valueFrom: |-
        ${
            if (inputs.memory_per_job) {
                return '\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\"';
            }
            return '\"-Xmx2048M\"';
        }
      shellQuote: false
    - position: 3
      valueFrom: SetNmMdAndUqTags
      shellQuote: false
    - prefix: ''
      position: 4
      valueFrom: |-
        ${
            var tmp = [].concat(inputs.in_alignments);
            var ext = ""; 
            if (inputs.output_file_format) {
                ext = inputs.output_file_format;
            } else {
                ext = tmp[0].path.split('.').pop();
            }
            
            if (inputs.output_prefix) {
                return '-O ' +  inputs.output_prefix + ".fixed." + ext;
            } else if (tmp[0].metadata && tmp[0].metadata.sample_id) {
                return '-O ' +  tmp[0].metadata.sample_id + ".fixed." + ext;
            } else {
                return '-O ' +  tmp[0].path.split('/').pop().split(".")[0] + ".fixed." + ext;
            }
            
        }
      shellQuote: false
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-setnmmdanduqtags-4-1-0-0/10
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    sbg:content_hash: a31d48359c8ea5e8ac91b2096488ac9e8a71d49dd3aa1a8ffbdcc09665a2c1f39
    sbg:contributors:
    - nens
    - uros_sipetic
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/15
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1555498307
    sbg:id: h-c05a5b7a/h-3d555f9f/h-25438681/0
    sbg:image_url:
    sbg:latestRevision: 10
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - id: https://software.broadinstitute.org/gatk/
      label: Homepage
    - id: https://github.com/broadinstitute/gatk/
      label: Source Code
    - id: |-
        https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
      label: Download
    - id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
      label: Publications
    - id: |-
        https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_SetNmMdAndUqTags.php
      label: Documentation
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1558518048
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 10
    sbg:revisionNotes: |-
      Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/15
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1555498307
      sbg:revision: 0
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/1
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1555582274
      sbg:revision: 1
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/5
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1556194603
      sbg:revision: 2
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/6
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557399646
      sbg:revision: 3
      sbg:revisionNotes: app info improved - perf bench needed
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557417063
      sbg:revision: 4
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/7
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557734531
      sbg:revision: 5
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/9
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000576
      sbg:revision: 6
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/10
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558100350
      sbg:revision: 7
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/11
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351574
      sbg:revision: 8
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/13
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558450064
      sbg:revision: 9
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/14
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558518048
      sbg:revision: 10
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/15
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: out_alignments
  sbg:x: 675.0732421875
  sbg:y: 260.1669006347656
- id: gatk_baserecalibrator_4_1_0_0
  label: GATK BaseRecalibrator
  in:
  - id: in_alignments
    source:
    - gatk_setnmmdanduqtags_4_1_0_0/out_alignments
  - id: include_intervals_file
    source: sbg_lines_to_interval_list_br/out_intervals
  - id: known_sites
    source:
    - known_sites
  - id: in_reference
    source: in_reference
  - id: use_original_qualities
    default: true
  scatter:
  - include_intervals_file
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK BaseRecalibrator CWL1.0
    doc: |-
      **GATK BaseRecalibrator** generates a recalibration table based on various covariates for input mapped read data [1]. It performs the first pass of the Base Quality Score Recalibration (BQSR) by assessing base quality scores of the input data.

      *A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*

      ###Common Use Cases

      * The **GATK BaseRecalibrator** tool requires the input mapped read data whose quality scores need to be assessed on its **Input alignments** (`--input`) input, the database of known polymorphic sites to skip over on its **Known sites** (`--known-sites`) input and a reference file on its **Reference** (`--reference`) input. On its **Output recalibration report** output, the tool generates a GATK report file with many tables: the list of arguments, the quantized qualities table, the recalibration table by read group, the recalibration table by quality score,
      the recalibration table for all the optional covariates [1].

      * Usage example:

      ```
      gatk --java-options "-Xmx2048M" BaseRecalibrator \
         --input my_reads.bam \
         --reference reference.fasta \
         --known-sites sites_of_variation.vcf \
         --known-sites another/optional/setOfSitesToMask.vcf \
         --output recal_data.table

      ```

      ###Changes Introduced by Seven Bridges

      * The output file will be prefixed using the **Output name prefix** parameter. If this value is not set, the output name will be generated based on the **Sample ID** metadata value from the input alignment file. If the **Sample ID** value is not set, the name will be inherited from the input alignment file name. In case there are multiple files on the **Input alignments** input, the files will be sorted by name and output file name will be generated based on the first file in the sorted file list, following the rules defined in the previous case. Moreover,  **recal_data** will be added before the extension of the output file name which is **CSV** by default.

      * **Include intervals** (`--intervals`) option is divided into **Include intervals string** and **Include intervals file** options.

      * **Exclude intervals** (`--exclude-intervals`) option is divided into **Exclude intervals string** and **Exclude intervals file** options.

      * The following GATK parameters were excluded from the tool wrapper: `--add-output-sam-program-record`, `--add-output-vcf-command-line`, `--arguments_file`, `--cloud-index-prefetch-buffer`, `--cloud-prefetch-buffer`, `--create-output-bam-index`, `--create-output-bam-md5`, `--create-output-variant-index`, `--create-output-variant-md5`, `--gatk-config-file`, `--gcs-max-retries`, `--gcs-project-for-requester-pays`, `--help`, `--lenient`, `--QUIET`, `--sites-only-vcf-output`, `--showHidden`, `--tmp-dir`, `--use-jdk-deflater`, `--use-jdk-inflater`, `--verbosity`, `--version`



      ###Common Issues and Important Notes

      *  **Memory per job** (`mem_per_job`) input allows a user to set the desired memory requirement when running a tool or adding it to a workflow. This input should be defined in MB. It is propagated to the Memory requirements part and “-Xmx” parameter of the tool. The default value is 2048MB.
      * **Memory overhead per job** (`mem_overhead_per_job`) input allows a user to set the desired overhead memory when running a tool or adding it to a workflow. This input should be defined in MB. This amount will be added to the Memory per job in the Memory requirements section but it will not be added to the “-Xmx” parameter. The default value is 100MB. 
      * Note: GATK tools that take in mapped read data expect a BAM file as the primary format [2]. More on GATK requirements for mapped sequence data formats can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats).
      * Note: **Known sites**, **Input alignments** should have corresponding index files in the same folder. 
      * Note: **Reference** FASTA file should have corresponding .fai (FASTA index) and .dict (FASTA dictionary) files in the same folder. 
      * Note: These **Read Filters** (`--read-filter`) are automatically applied to the data by the Engine before processing by **BaseRecalibrator** [1]: **NotSecondaryAlignmentReadFilter**, **PassesVendorQualityCheckReadFilter**, **MappedReadFilter**, **MappingQualityAvailableReadFilter**, **NotDuplicateReadFilter**, **MappingQualityNotZeroReadFilter**, **WellformedReadFilter**
      * Note: If the **Read filter** (`--read-filter`) option is set to "LibraryReadFilter", the **Library** (`--library`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "PlatformReadFilter", the **Platform filter name** (`--platform-filter-name`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to"PlatformUnitReadFilter", the **Black listed lanes** (`--black-listed-lanes`) option must be set to some value. 
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadGroupBlackListReadFilter", the **Read group black list** (`--read-group-black-list`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadGroupReadFilter", the **Keep read group** (`--keep-read-group`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadLengthReadFilter", the **Max read length** (`--max-read-length`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadNameReadFilter", the **Read name** (`--read-name`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadStrandFilter", the **Keep reverse strand only** (`--keep-reverse-strand-only`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "SampleReadFilter", the **Sample** (`--sample`) option must be set to some value.
      * Note: The following options are valid only if the appropriate **Read filter** (`--read-filter`) is specified: **Ambig filter bases** (`--ambig-filter-bases`), **Ambig filter frac** (`--ambig-filter-frac`), **Max fragment length** (`--max-fragment-length`), **Maximum mapping quality** (`--maximum-mapping-quality`), **Minimum mapping quality** (`--minimum-mapping-quality`),  **Do not require soft clips** (`--dont-require-soft-clips-both-ends`), **Filter too short** (`--filter-too-short`), **Min read length** (`--min-read-length`). See the description of each parameter for information on the associated **Read filter**.
      * Note: The wrapper has not been tested for the SAM file type on the **Input alignments** input port, nor for the BCF file type on the **Known sites** input port.

      ###Performance Benchmarking

      Below is a table describing runtimes and task costs of **GATK BaseRecalibrator** for a couple of different samples, executed on AWS cloud instances:

      | Experiment type |  Input size | Duration |  Cost (on-demand) | Instance (AWS) | 
      |:--------------:|:------------:|:--------:|:-------:|:---------:|
      |     RNA-Seq     |  2.2 GB |   9min   | ~0.08$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     |  6.6 GB |   19min   | ~0.17$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 11 GB |  27min  | ~0.24$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 22 GB |  46min  | ~0.41$ | c4.2xlarge (8 CPUs) |

      *Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*

      ###References

      [1] [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036726891-BaseRecalibrator)

      [2] [GATK Mapped sequence data formats](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: "${\n    return inputs.cpu_per_job ? inputs.cpu_per_job : 1;\n}"
      ramMin: |-
        ${
          var memory = 2048;
          
          if(inputs.mem_per_job) {
          	 memory = inputs.mem_per_job;
          }
          if(inputs.mem_overhead_per_job) {
        	memory += inputs.mem_overhead_per_job;
          }
          else {
             memory += 100;
          }
          return memory;
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/marijeta_slavkovic/gatk-4-1-0-0:0
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib:
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file)) {
                file['metadata'] = {}
            }
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
            return file
        };
        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!o2) {
                return o1;
            };
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
                for (var key in commonMetadata) {
                    if (!(key in example)) {
                        delete commonMetadata[key]
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
                if (o1.secondaryFiles) {
                    o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)
                }
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                    if (o1[i].secondaryFiles) {
                        o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)
                    }
                }
            }
            return o1;
        };

    inputs:
    - id: ambig_filter_bases
      label: Ambig filter bases
      doc: |-
        Valid only if "AmbiguousBaseReadFilter" is specified:
        Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction. Cannot be used in conjuction with argument(s) ambig-filter-frac.
      type: int?
      inputBinding:
        prefix: --ambig-filter-bases
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: 'null'
    - id: ambig_filter_frac
      label: Ambig filter frac
      doc: |-
        Valid only if "AmbiguousBaseReadFilter" is specified:
        Threshold fraction of ambiguous bases. Cannot be used in conjuction with argument(s) ambig-filter-bases.
      type: float?
      inputBinding:
        prefix: --ambig-filter-frac
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '0.05'
    - id: binary_tag_name
      label: Binary tag name
      doc: The binary tag covariate name if using it.
      type: string?
      inputBinding:
        prefix: --binary-tag-name
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: black_listed_lanes
      label: Black listed lanes
      doc: |-
        Valid only if "PlatformUnitReadFilter" is specified:
        Platform unit (PU) to filter out. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.black_listed_lanes)
              {
                  var bl_lanes = [].concat(inputs.black_listed_lanes);
                  var cmd = [];
                  for (var i = 0; i < bl_lanes.length; i++) 
                  {
                      cmd.push('--black-listed-lanes', bl_lanes[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: bqsr_baq_gap_open_penalty
      label: BQSR BAQ gap open penalty
      doc: |-
        BQSR BAQ gap open penalty (Phred Scaled). Default value is 40. 30 is perhaps better for whole genome call sets.
      type: float?
      inputBinding:
        prefix: --bqsr-baq-gap-open-penalty
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '40'
    - id: default_base_qualities
      label: Default base qualities
      doc: Assign a default base quality.
      type: int?
      inputBinding:
        prefix: --default-base-qualities
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '-1'
    - id: deletions_default_quality
      label: Deletions default quality
      doc: Default quality for the base deletions covariate.
      type: int?
      inputBinding:
        prefix: --deletions-default-quality
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '45'
    - id: disable_read_filter
      label: Disable read filter
      doc: |-
        Read filters to be disabled before analysis. This argument may be specified 0 or more times.
      type:
      - 'null'
      - type: array
        items:
          name: disable_read_filter
          type: enum
          symbols:
          - MappedReadFilter
          - MappingQualityAvailableReadFilter
          - MappingQualityNotZeroReadFilter
          - NotDuplicateReadFilter
          - NotSecondaryAlignmentReadFilter
          - PassesVendorQualityCheckReadFilter
          - WellformedReadFilter
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--disable-read-filter', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:altPrefix: -DF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: disable_sequence_dictionary_validation
      label: Disable sequence dictionary validation
      doc: |-
        If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!
      type: boolean?
      inputBinding:
        prefix: --disable-sequence-dictionary-validation
        position: 4
        shellQuote: false
      sbg:altPrefix: -disable-sequence-dictionary-validation
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: disable_tool_default_read_filters
      label: Disable tool default read filters
      doc: |-
        Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on).
      type: boolean?
      inputBinding:
        prefix: --disable-tool-default-read-filters
        position: 4
        shellQuote: false
      sbg:altPrefix: -disable-tool-default-read-filters
      sbg:category: Advanced Arguments
      sbg:toolDefaultValue: 'false'
    - id: dont_require_soft_clips_both_ends
      label: Dont require soft clips both ends
      doc: |-
        Valid only if "OverclippedReadFilter" is specified:
        Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block.
      type: boolean?
      inputBinding:
        prefix: --dont-require-soft-clips-both-ends
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: 'false'
    - id: exclude_intervals_file
      label: Exclude intervals file
      doc: One or more genomic intervals to exclude from processing.
      type: File?
      inputBinding:
        prefix: --exclude-intervals
        position: 4
        shellQuote: false
      sbg:altPrefix: -XL
      sbg:category: Optional Arguments
      sbg:fileTypes: BED, LIST, INTERVAL_LIST
      sbg:toolDefaultValue: 'null'
    - id: exclude_intervals_string
      label: Exclude intervals string
      doc: |-
        One or more genomic intervals to exclude from processing. This argument may be specified 0 or more times.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |+
          ${
              if (inputs.exclude_intervals_string)
              {
                  var exclude_string = [].concat(inputs.exclude_intervals_string);
                  var cmd = [];
                  for (var i = 0; i < exclude_string.length; i++) 
                  {
                      cmd.push('--exclude-intervals', exclude_string[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }


        shellQuote: false
      sbg:altPrefix: -XL
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: filter_too_short
      label: Filter too short
      doc: |-
        Valid only if "OverclippedReadFilter" is specified:
        Minimum number of aligned bases.
      type: int?
      inputBinding:
        prefix: --filter-too-short
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '30'
    - id: indels_context_size
      label: Indels context size
      doc: Size of the k-mer context to be used for base insertions and deletions.
      type: int?
      inputBinding:
        prefix: --indels-context-size
        position: 4
        shellQuote: false
      sbg:altPrefix: -ics
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '3'
    - id: in_alignments
      label: Input alignments
      doc: |-
        BAM/SAM/CRAM file containing reads. This argument must be specified at least once.
      type: File[]
      secondaryFiles:
      - |-
        ${
            var in_alignments = self;
            if (in_alignments.nameext == '.bam' || in_alignments.nameext == '.BAM') {
                return [in_alignments.basename + ".bai", in_alignments.nameroot + ".bai"];
            }
            else if (in_alignments.nameext == '.cram' || in_alignments.nameext == '.CRAM') {
                return [in_alignments.basename + ".crai", in_alignments.nameroot + ".crai", in_alignments.basename + ".bai"];     
            }
            return '';
        }
      inputBinding:
        position: 4
        valueFrom: |
          ${
              if (inputs.in_alignments) {
                  var alignments = [].concat(inputs.in_alignments);
                  var cmd = [];
                  for (var i=0; i<alignments.length; i++) {
                      cmd.push('--input', alignments[i].path);
                  }
                  return cmd.join(' ');
              } 
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, CRAM
    - id: insertions_default_quality
      label: Insertions default quality
      doc: Default quality for the base insertions covariate.
      type: int?
      inputBinding:
        prefix: --insertions-default-quality
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '45'
    - id: interval_exclusion_padding
      label: Interval exclusion padding
      doc: Amount of padding (in bp) to add to each interval you are excluding.
      type: int?
      inputBinding:
        prefix: --interval-exclusion-padding
        position: 4
        shellQuote: false
      sbg:altPrefix: -ixp
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: interval_merging_rule
      label: Interval merging rule
      doc: Interval merging rule for abutting intervals.
      type:
      - 'null'
      - name: interval_merging_rule
        type: enum
        symbols:
        - ALL
        - OVERLAPPING_ONLY
      inputBinding:
        prefix: --interval-merging-rule
        position: 4
        shellQuote: false
      sbg:altPrefix: -imr
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: ALL
    - id: interval_padding
      label: Interval padding
      doc: Amount of padding (in bp) to add to each interval you are including.
      type: int?
      inputBinding:
        prefix: --interval-padding
        position: 4
        shellQuote: false
      sbg:altPrefix: -ip
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: interval_set_rule
      label: Interval set rule
      doc: Set merging approach to use for combining interval inputs.
      type:
      - 'null'
      - name: interval_set_rule
        type: enum
        symbols:
        - UNION
        - INTERSECTION
      inputBinding:
        prefix: --interval-set-rule
        position: 4
        shellQuote: false
      sbg:altPrefix: -isr
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: UNION
    - id: include_intervals_file
      label: Include intervals file
      doc: One or more genomic intervals over which to operate.
      type: File?
      inputBinding:
        prefix: --intervals
        position: 4
        shellQuote: false
      sbg:altPrefix: -L
      sbg:category: Optional Arguments
      sbg:fileTypes: BED, LIST, INTERVAL_LIST
      sbg:toolDefaultValue: 'null'
    - id: include_intervals_string
      label: Include intervals string
      doc: |-
        One or more genomic intervals over which to operate. This argument may be specified 0 or more times.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |+
          ${
              if (inputs.include_intervals_string)
              {
                  var include_string = [].concat(inputs.include_intervals_string);
                  var cmd = [];
                  for (var i = 0; i < include_string.length; i++) 
                  {
                      cmd.push('--intervals', include_string[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }


        shellQuote: false
      sbg:altPrefix: -L
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: keep_read_group
      label: Keep read group
      doc: |-
        Valid only if "ReadGroupReadFilter" is specified:
        The name of the read group to keep. Required.
      type: string?
      inputBinding:
        prefix: --keep-read-group
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: keep_reverse_strand_only
      label: Keep reverse strand only
      doc: |-
        Valid only if "ReadStrandFilter" is specified:
        Keep only reads on the reverse strand. Required.
      type:
      - 'null'
      - name: keep_reverse_strand_only
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --keep-reverse-strand-only
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: known_sites
      label: Known sites
      doc: |-
        One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.  This argument must be specified at least once.
      type: File[]
      secondaryFiles:
      - |-
        ${
            var in_sites = self;
            if (in_sites.nameext == ".gz" || in_sites.nameext == '.GZ') {
                    var tmp = in_sites.basename.slice(-7);
                    if(tmp.toLowerCase() == '.vcf.gz') {
                        return in_sites.basename + ".tbi";  
                    }
            }
            else if (in_sites.nameext == '.vcf' || in_sites.nameext == '.VCF' || in_sites.nameext == '.bed' || in_sites.nameext == '.BED') {
                return in_sites.basename + ".idx";
            }
            return in_sites.basename + ".idx";
        }
      inputBinding:
        position: 5
        valueFrom: |-
          ${
              if (inputs.known_sites)
              {
                  var sites = [].concat(inputs.known_sites);
                  var cmd = [];
                  for (var i = 0; i < sites.length; i++) 
                  {
                      cmd.push('--known-sites', sites[i].path);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:category: Required Arguments
      sbg:fileTypes: VCF, VCF.GZ, BED
    - id: library
      label: Library
      doc: |-
        Valid only if "LibraryReadFilter" is specified:
        Name of the library to keep. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.library)
              {
                  var lib = [].concat(inputs.library);
                  var cmd = [];
                  for (var i = 0; i < lib.length; i++) 
                  {
                      cmd.push('--library', lib[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -library
      sbg:category: Conditional Arguments for readFilter
    - id: low_quality_tail
      label: Low quality tail
      doc: Minimum quality for the bases in the tail of the reads to be considered.
      type: int?
      inputBinding:
        prefix: --low-quality-tail
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
    - id: max_fragment_length
      label: Max fragment length
      doc: |-
        Valid only if "FragmentLengthReadFilter" is specified:
        Maximum length of fragment (insert size).
      type: int?
      inputBinding:
        prefix: --max-fragment-length
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '1000000'
    - id: max_read_length
      label: Max read length
      doc: |-
        Valid only if "ReadLengthReadFilter" is specified:
        Keep only reads with length at most equal to the specified value. Required.
      type: int?
      inputBinding:
        prefix: --max-read-length
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: maximum_cycle_value
      label: Maximum cycle value
      doc: The maximum cycle value permitted for the Cycle covariate.
      type: int?
      inputBinding:
        prefix: --maximum-cycle-value
        position: 4
        shellQuote: false
      sbg:altPrefix: -max-cycle
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500'
    - id: maximum_mapping_quality
      label: Maximum mapping quality
      doc: |-
        Valid only if "MappingQualityReadFilter" is specified:
        Maximum mapping quality to keep (inclusive).
      type: int?
      inputBinding:
        prefix: --maximum-mapping-quality
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: 'null'
    - id: mem_overhead_per_job
      label: Memory overhead per job
      doc: |-
        It allows a user to set the desired overhead memory (in MB) when running a tool or adding it to a workflow.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '100'
    - id: mem_per_job
      label: Memory per job
      doc: |-
        It allows a user to set the desired memory requirement (in MB) when running a tool or adding it to a workflow.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '2048'
    - id: min_read_length
      label: Min read length
      doc: |-
        Valid only if "ReadLengthReadFilter" is specified:
        Keep only reads with length at least equal to the specified value.
      type: int?
      inputBinding:
        prefix: --min-read-length
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '1'
    - id: minimum_mapping_quality
      label: Minimum mapping quality
      doc: |-
        Valid only if "MappingQualityReadFilter" is specified:
        Minimum mapping quality to keep (inclusive).
      type: int?
      inputBinding:
        prefix: --minimum-mapping-quality
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '10'
    - id: mismatches_context_size
      label: Mismatches context size
      doc: Size of the k-mer context to be used for base mismatches.
      type: int?
      inputBinding:
        prefix: --mismatches-context-size
        position: 4
        shellQuote: false
      sbg:altPrefix: -mcs
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
    - id: mismatches_default_quality
      label: Mismatches default quality
      doc: Default quality for the base mismatches covariate.
      type: int?
      inputBinding:
        prefix: --mismatches-default-quality
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '-1'
    - id: platform_filter_name
      label: Platform filter name
      doc: |-
        Valid only if "PlatformReadFilter" is specified:
        Platform attribute (PL) to match. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.platform_filter_name)
              {
                  var pfn = [].concat(inputs.platform_filter_name);
                  var cmd = [];
                  for (var i = 0; i < pfn.length; i++) 
                  {
                      cmd.push('--platform-filter-name', pfn[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: preserve_qscores_less_than
      label: Preserve qscores less than
      doc: |-
        Don't recalibrate bases with quality scores less than this threshold (with -bqsr).
      type: int?
      inputBinding:
        prefix: --preserve-qscores-less-than
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '6'
    - id: quantizing_levels
      label: Quantizing levels
      doc: Number of distinct quality scores in the quantized output.
      type: int?
      inputBinding:
        prefix: --quantizing-levels
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '16'
    - id: read_filter
      label: Read filter
      doc: |-
        Read filters to be applied before analysis. This argument may be specified 0 or more times.
      type:
      - 'null'
      - type: array
        items:
          name: read_filter
          type: enum
          symbols:
          - AlignmentAgreesWithHeaderReadFilter
          - AllowAllReadsReadFilter
          - AmbiguousBaseReadFilter
          - CigarContainsNoNOperator
          - FirstOfPairReadFilter
          - FragmentLengthReadFilter
          - GoodCigarReadFilter
          - HasReadGroupReadFilter
          - LibraryReadFilter
          - MappedReadFilter
          - MappingQualityAvailableReadFilter
          - MappingQualityNotZeroReadFilter
          - MappingQualityReadFilter
          - MatchingBasesAndQualsReadFilter
          - MateDifferentStrandReadFilter
          - MateOnSameContigOrNoMappedMateReadFilter
          - MetricsReadFilter
          - NonChimericOriginalAlignmentReadFilter
          - NonZeroFragmentLengthReadFilter
          - NonZeroReferenceLengthAlignmentReadFilter
          - NotDuplicateReadFilter
          - NotOpticalDuplicateReadFilter
          - NotSecondaryAlignmentReadFilter
          - NotSupplementaryAlignmentReadFilter
          - OverclippedReadFilter
          - PairedReadFilter
          - PassesVendorQualityCheckReadFilter
          - PlatformReadFilter
          - PlatformUnitReadFilter
          - PrimaryLineReadFilter
          - ProperlyPairedReadFilter
          - ReadGroupBlackListReadFilter
          - ReadGroupReadFilter
          - ReadLengthEqualsCigarLengthReadFilter
          - ReadLengthReadFilter
          - ReadNameReadFilter
          - ReadStrandFilter
          - SampleReadFilter
          - SecondOfPairReadFilter
          - SeqIsStoredReadFilter
          - ValidAlignmentEndReadFilter
          - ValidAlignmentStartReadFilter
          - WellformedReadFilter
      inputBinding:
        prefix: ''
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--read-filter', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:altPrefix: -RF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: read_group_black_list
      label: Read group black list
      doc: |-
        Valid only if "ReadGroupBlackListReadFilter" is specified:
        The name of the read group to filter out. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.read_group_black_list)
              {
                  var rgbl = [].concat(inputs.read_group_black_list);
                  var cmd = [];
                  for (var i = 0; i < rgbl.length; i++) 
                  {
                      cmd.push('--read-group-black-list', rgbl[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: read_name
      label: Read name
      doc: |-
        Valid only if "ReadNameReadFilter" is specified:
        Keep only reads with this read name. Required.
      type: string?
      inputBinding:
        prefix: --read-name
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: read_validation_stringency
      label: Read validation stringency
      doc: |-
        Validation stringency for all SAM/BAM/CRAM/SRA files read by this program. The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
      type:
      - 'null'
      - name: read_validation_stringency
        type: enum
        symbols:
        - STRICT
        - LENIENT
        - SILENT
      inputBinding:
        prefix: --read-validation-stringency
        position: 4
        shellQuote: false
      sbg:altPrefix: -VS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: SILENT
    - id: in_reference
      label: Reference
      doc: Reference sequence file.
      type: File
      secondaryFiles:
      - .fai
      - ^.dict
      inputBinding:
        prefix: --reference
        position: 4
        shellQuote: false
      sbg:altPrefix: -R
      sbg:category: Required Arguments
      sbg:fileTypes: FASTA, FA
    - id: sample
      label: Sample
      doc: |-
        Valid only if "SampleReadFilter" is specified:
        The name of the sample(s) to keep, filtering out all others. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.sample)
              {
                  var samp = [].concat(inputs.sample);
                  var cmd = [];
                  for (var i = 0; i < samp.length; i++) 
                  {
                      cmd.push('--sample', samp[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -sample
      sbg:category: Conditional Arguments for readFilter
    - id: sequence_dictionary
      label: Sequence dictionary
      doc: |-
        Use the given sequence dictionary as the master/canonical sequence dictionary. Must be a .dict file.
      type: File?
      inputBinding:
        prefix: --sequence-dictionary
        position: 4
        shellQuote: false
      sbg:altPrefix: -sequence-dictionary
      sbg:category: Optional Arguments
      sbg:fileTypes: DICT
      sbg:toolDefaultValue: '10.0'
    - id: use_original_qualities
      label: Use original qualities
      doc: Use the base quality scores from the OQ tag.
      type: boolean?
      inputBinding:
        prefix: --use-original-qualities
        position: 4
        shellQuote: false
      sbg:altPrefix: -OQ
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: prefix
      label: Output name prefix
      doc: Output file name prefix.
      type: string?
      sbg:category: Config Inputs
    - id: cpu_per_job
      label: CPU per job
      doc: CPU per job.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
    - id: disable_bam_index_caching
      label: Disable BAM index caching
      doc: |-
        If true, don't cache BAM indexes, this will reduce memory requirements but may harm performance if many intervals are specified. Caching is automatically disabled if there are no intervals specified.
      type: boolean?
      inputBinding:
        prefix: --disable-bam-index-caching
        position: 4
        shellQuote: false
      sbg:altPrefix: -DBIC
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: seconds_between_progress_updates
      label: Seconds between progress updates
      doc: Output traversal statistics every time this many seconds elapse.
      type: float?
      inputBinding:
        prefix: --seconds-between-progress-updates
        position: 4
        shellQuote: false
      sbg:altPrefix: -seconds-between-progress-updates
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '10.00'
    - id: read_index
      label: Read index
      doc: |-
        Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically. This argument may be specified 0 or more times.
      type: File[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.read_index)
              {
                  var r_index = [].concat(inputs.read_index);
                  var cmd = [];
                  for (var i = 0; i < r_index.length; i++) 
                  {
                      cmd.push('--read-index', r_index[i].path);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -read-index
      sbg:category: Optional Arguments
      sbg:fileTypes: BAI, CRAI

    outputs:
    - id: out_bqsr_report
      label: Output recalibration report
      doc: The output recalibration table file to create.
      type: File?
      outputBinding:
        glob: '*.csv'
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: CSV

    baseCommand:
    - /opt/gatk-4.1.0.0/gatk --java-options
    arguments:
    - prefix: ''
      position: 1
      valueFrom: |-
        ${
            if (inputs.mem_per_job) {
                return '\"-Xmx'.concat(inputs.mem_per_job, 'M') + '\"';
            } else {
                return '\"-Xmx2048M\"';
            }
        }
      shellQuote: false
    - prefix: ''
      position: 2
      valueFrom: BaseRecalibrator
      shellQuote: false
    - prefix: --output
      position: 3
      valueFrom: |-
        ${
            //sort list of input files by nameroot
            function sortNameroot(x, y) {
                if (x.nameroot < y.nameroot) {
                    return -1;
                }
                if (x.nameroot > y.nameroot) {
                    return 1;
                }
                return 0;
            }
                
            var output_prefix;
            var in_num = [].concat(inputs.in_alignments).length;
            var in_align = [].concat(inputs.in_alignments);
            
            //if input_prefix is provided by the user
            if (inputs.prefix) {
                output_prefix = inputs.prefix;
                if (in_num > 1) {
                    output_prefix = output_prefix + '.' + in_num;
                }
            }
            else {
                //if there is only one input file
                if(in_num == 1){
                    // check if the sample_id metadata value is defined for the input file
                    if(in_align[0].metadata && in_align[0].metadata.sample_id) {
                        output_prefix = in_align[0].metadata.sample_id;
                    // if sample_id is not defined
                    } else {
                        output_prefix = in_align[0].path.split('/').pop().split('.')[0];
                    }
                }
                //if there are more than 1 input files
                //sort list of input file objects alphabetically by file name 
                //take the first element from that list, and generate output file name as if that file is the only file on the input. 
                else if(in_num > 1) {
                    //sort list of input files by nameroot
                    in_align.sort(sortNameroot);
                    //take the first alphabetically sorted file
                    var first_file = in_align[0];
                    //check if the sample_id metadata value is defined for the input file
                    if(first_file.metadata && first_file.metadata.sample_id) {
                        output_prefix = first_file.metadata.sample_id + '.' + in_num;
                    // if sample_id is not defined
                    } else {
                        output_prefix = first_file.path.split('/').pop().split('.')[0] + '.' + in_num;
                    }
                }
            }
            var output_full = output_prefix + '.recal_data.csv';
            return output_full;
        }
      shellQuote: false
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-baserecalibrator-4-1-0-0/22
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    - CWL1.0
    sbg:content_hash: af89c0ecbd011d6f1e94510e1c0947c9cce2b6d5d05713be641ff8cbc7de1d6af
    sbg:contributors:
    - nens
    - veliborka_josipovic
    - uros_sipetic
    - marijeta_slavkovic
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552922094
    sbg:id: h-9342c502/h-42dccb3d/h-ed8e8a69/0
    sbg:image_url:
    sbg:latestRevision: 22
    sbg:license: BSD 3-Clause License
    sbg:links:
    - id: https://www.broadinstitute.org/gatk/index.php
      label: Homepage
    - id: https://github.com/broadinstitute/gatk
      label: Source Code
    - id: |-
        https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
      label: Download
    - id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
      label: Publication
    - id: https://gatk.broadinstitute.org/hc/en-us/articles/360036726891-BaseRecalibrator
      label: Documentation
    sbg:modifiedBy: marijeta_slavkovic
    sbg:modifiedOn: 1603296363
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 22
    sbg:revisionNotes: |-
      secondary files known_sites (return basename.idx instead of '' when not VCF or VCF.GZ), small description
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1552922094
      sbg:revision: 0
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/11
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554492924
      sbg:revision: 1
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/14
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554492998
      sbg:revision: 2
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/15
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554720866
      sbg:revision: 3
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/17
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554999207
      sbg:revision: 4
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/18
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1556030757
      sbg:revision: 5
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/19
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557735256
      sbg:revision: 6
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/20
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000594
      sbg:revision: 7
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/21
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351546
      sbg:revision: 8
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/23
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558450805
      sbg:revision: 9
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/24
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558517350
      sbg:revision: 10
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/25
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558518057
      sbg:revision: 11
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/26
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1571321280
      sbg:revision: 12
      sbg:revisionNotes: known_snps null handled
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1593698771
      sbg:revision: 13
      sbg:revisionNotes: New wrapper
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1593699523
      sbg:revision: 14
      sbg:revisionNotes: Description review suggestions added
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1593699583
      sbg:revision: 15
      sbg:revisionNotes: Description review suggestions added
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1594047999
      sbg:revision: 16
      sbg:revisionNotes: naming description and benchmarking price review
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1594725435
      sbg:revision: 17
      sbg:revisionNotes: added CRAM and SAM to suggested types for in_alignments
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1594725563
      sbg:revision: 18
      sbg:revisionNotes: removed SAM as file suggestion
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1597669945
      sbg:revision: 19
      sbg:revisionNotes: changed default mem_per_job to 2048
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1598131454
      sbg:revision: 20
      sbg:revisionNotes: added [].concat to arrays
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1603199349
      sbg:revision: 21
      sbg:revisionNotes: description edited (usage example Xmx, memory in description
        etc)
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1603296363
      sbg:revision: 22
      sbg:revisionNotes: |-
        secondary files known_sites (return basename.idx instead of '' when not VCF or VCF.GZ), small description
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: out_bqsr_report
  sbg:x: 1241.2686767578125
  sbg:y: 307.5648193359375
- id: gatk_createsequencegroupingtsv_4_1_0_0
  label: GATK CreateSequenceGroupingTSV
  in:
  - id: ref_dict
    source: ref_dict
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK CreateSequenceGroupingTSV
    doc: |-
      **CreateSequenceGroupingTSV** tool generate sets of intervals for scatter-gathering over chromosomes.

      It takes **Reference dictionary** file (`--ref_dict`) as an input and creates files which contain chromosome names grouped based on their sizes.


      ###**Common Use Cases**

      The tool has only one input (`--ref_dict`) which is required and has no additional arguments. **CreateSequenceGroupingTSV** tool results are **Sequence Grouping** file which is a text file containing chromosome groups, and **Sequence Grouping with Unmapped**, a text file which has the same content as **Sequence Grouping** with additional line containing "unmapped" string.


      * Usage example


      ```
      python CreateSequenceGroupingTSV.py 
            --ref_dict example_reference.dict

      ```



      ###**Changes Introduced by Seven Bridges**

      Python code provided within WGS Germline WDL was adjusted to be called as a script (`CreateSequenceGroupingTSV.py`).


      ###**Common Issues and Important Notes**

      None.


      ### Reference
      [1] [CreateSequenceGroupingTSV](https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 1000
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0
    - class: InitialWorkDirRequirement
      listing:
      - entryname: CreateSequenceGroupingTSV.py
        writable: false
        entry: |-
          import argparse

          args = argparse.ArgumentParser(description='This tool takes reference dictionary file as an input'
                                                       ' and creates files which contain chromosome names grouped'
                                                       ' based on their sizes.')

          args.add_argument('--ref_dict', help='Reference dictionary', required=True)
          parsed = args.parse_args()
          ref_dict = parsed.ref_dict

          with open(ref_dict, 'r') as ref_dict_file:
              sequence_tuple_list = []
              longest_sequence = 0
              for line in ref_dict_file:
                  if line.startswith("@SQ"):
                      line_split = line.split("\t")
                      # (Sequence_Name, Sequence_Length)
                      sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
              longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
          # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
          # the last element after a :, so we add this as a sacrificial element.
          hg38_protection_tag = ":1+"
          # initialize the tsv string with the first sequence
          tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
          temp_size = sequence_tuple_list[0][1]
          for sequence_tuple in sequence_tuple_list[1:]:
              if temp_size + sequence_tuple[1] <= longest_sequence:
                  temp_size += sequence_tuple[1]
                  tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
              else:
                  tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
                  temp_size = sequence_tuple[1]
          # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
          with open("./sequence_grouping.txt", "w") as tsv_file:
              tsv_file.write(tsv_string)
              tsv_file.close()

          tsv_string += '\n' + "unmapped"

          with open("./sequence_grouping_with_unmapped.txt", "w") as tsv_file_with_unmapped:
              tsv_file_with_unmapped.write(tsv_string)
              tsv_file_with_unmapped.close()
    - class: InlineJavascriptRequirement
      expressionLib:
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

    inputs:
    - id: ref_dict
      label: Reference Dictionary
      doc: |-
        Reference dictionary containing information about chromosome names and their lengths.
      type: File
      inputBinding:
        prefix: --ref_dict
        position: 0
        shellQuote: false
      sbg:fileTypes: DICT

    outputs:
    - id: sequence_grouping
      label: Sequence Grouping
      doc: |-
        Each line of the file represents one group of chromosomes which are processed together in later steps of the GATK Germline workflow. The groups are determined based on the chromosomes sizes.
      type: File?
      outputBinding:
        glob: sequence_grouping.txt
        outputEval: $(inheritMetadata(self, inputs.ref_dict))
      sbg:fileTypes: TXT
    - id: sequence_grouping_with_unmapped
      label: Sequence Grouping with Unmapped
      doc: |-
        The file has the same content as "Sequence Grouping" file, with an additional, last line containing "unmapped" string.
      type: File?
      outputBinding:
        glob: sequence_grouping_with_unmapped.txt
        outputEval: $(inheritMetadata(self, inputs.ref_dict))
      sbg:fileTypes: TXT

    baseCommand:
    - python
    - CreateSequenceGroupingTSV.py
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-createsequencegroupingtsv-4-1-0-0/4
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BED Processing
    sbg:content_hash: a9afa170a339934c60906ff616a6f2155426a9df80067bfc64f4140593aeffda6
    sbg:contributors:
    - nens
    - uros_sipetic
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/6
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1555580154
    sbg:id: h-0e394414/h-f3dcf09c/h-c405fcdf/0
    sbg:image_url:
    sbg:latestRevision: 4
    sbg:license: BSD 3-clause
    sbg:links:
    - id: https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels
      label: GATK Germline GitHub
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1558351560
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 4
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/6
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1555580154
      sbg:revision: 0
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/1
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557734537
      sbg:revision: 1
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/3
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557914517
      sbg:revision: 2
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/4
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000609
      sbg:revision: 3
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/5
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351560
      sbg:revision: 4
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/6
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: sequence_grouping
  - id: sequence_grouping_with_unmapped
  sbg:x: 767.7706909179688
  sbg:y: 6.801900386810303
- id: gatk_gatherbqsrreports_4_1_0_0
  label: GATK GatherBQSRReports
  in:
  - id: in_bqsr_reports
    source:
    - gatk_baserecalibrator_4_1_0_0/out_bqsr_report
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK GatherBQSRReports CWL1.0
    doc: |-
      **GATK GatherBQSRReports** gathers scattered BQSR recalibration reports into a single file [1].

      *A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*


      ### Common Use Cases 

      * This tool is intended to be used to combine recalibration tables from runs of **GATK BaseRecalibrator** parallelized per-interval.

      * Usage example:
      ```
         gatk --java-options "-Xmx2048M" GatherBQSRReports \
         --input input1.csv \
         --input input2.csv \
         --output output.csv

      ```


      ###Changes Introduced by Seven Bridges

      * The output file will be prefixed using the **Output name prefix** parameter. If this value is not set, the output name will be generated based on the **Sample ID** metadata value from **Input BQSR reports**. If the **Sample ID** value is not set, the name will be inherited from the **Input BQSR reports** file name. In case there are multiple files on the **Input BQSR reports** input, the files will be sorted by name and output file name will be generated based on the first file in the sorted file list, following the rules defined in the previous case. Moreover, **.recal_data** will be added before the extension of the output file name.

      * The following GATK parameters were excluded from the tool wrapper: `--arguments_file`, `--gatk-config-file`, `--gcs-max-retries`, `--gcs-project-for-requester-pays`, `--help`, `--QUIET`, `--showHidden`, `--tmp-dir`, `--use-jdk-deflater`, `--use-jdk-inflater`, `--verbosity`, `--version`


      ###Common Issues and Important Notes

      *  **Memory per job** (`mem_per_job`) input allows a user to set the desired memory requirement when running a tool or adding it to a workflow. This input should be defined in MB. It is propagated to the Memory requirements part and “-Xmx” parameter of the tool. The default value is 2048MB.

      * **Memory overhead per job** (`mem_overhead_per_job`) input allows a user to set the desired overhead memory when running a tool or adding it to a workflow. This input should be defined in MB. This amount will be added to the Memory per job in the Memory requirements section but it will not be added to the “-Xmx” parameter. The default value is 100MB. 


      ###Performance Benchmarking

      This tool is fast, with a running time of a few minutes. The experiment task was performed on the default AWS on-demand c4.2xlarge instance on 50 CSV files (size of each ~350KB) and took 2 minutes to finish ($0.02).

      *Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*


      ###References

      [1] [GATK GatherBQSRReports](https://gatk.broadinstitute.org/hc/en-us/articles/360036359192-GatherBQSRReports)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: '$(inputs.cpu_per_job ? inputs.cpu_per_job : 1)'
      ramMin: |-
        ${
          var memory = 2048;
          
          if(inputs.mem_per_job) {
          	 memory = inputs.mem_per_job;
          }
          if(inputs.mem_overhead_per_job) {
        	memory += inputs.mem_overhead_per_job;
          }
          else {
              memory += 100;
          }
          return memory;
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/marijeta_slavkovic/gatk-4-1-0-0:0
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib:
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file)) {
                file['metadata'] = {}
            }
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
            return file
        };
        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!o2) {
                return o1;
            };
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
                for (var key in commonMetadata) {
                    if (!(key in example)) {
                        delete commonMetadata[key]
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
                if (o1.secondaryFiles) {
                    o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)
                }
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                    if (o1[i].secondaryFiles) {
                        o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)
                    }
                }
            }
            return o1;
        };

    inputs:
    - id: in_bqsr_reports
      label: Input BQSR reports
      doc: |-
        List of scattered BQSR report files. This argument must be specified at least once.
      type: File[]
      inputBinding:
        position: 4
        valueFrom: |-
          ${
             if (inputs.in_bqsr_reports)
             {
                 var bqsr_reports = [].concat(inputs.in_bqsr_reports);
                 var cmd = [];
                 for (var i = 0; i < bqsr_reports.length; i++)
                 {
                     cmd.push('--input', bqsr_reports[i].path);
                 }
                 return cmd.join(' ');
             }
             return '';
          }
        itemSeparator: 'null'
        shellQuote: false
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: CSV
    - id: prefix
      label: Output name prefix
      doc: Output prefix for the gathered BQSR report.
      type: string?
      sbg:category: Config Inputs
    - id: mem_per_job
      label: Memory per job
      doc: |-
        It allows a user to set the desired memory requirement (in MB) when running a tool or adding it to a workflow.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '2048'
    - id: mem_overhead_per_job
      label: Memory overhead per job
      doc: |-
        It allows a user to set the desired overhead memory (in MB) when running a tool or adding it to a workflow.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '100'
    - id: cpu_per_job
      label: CPU per job
      doc: Number of CPUs to be used per job.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'

    outputs:
    - id: out_gathered_bqsr_report
      label: Gathered BQSR report
      doc: File to output the gathered file to.
      type: File?
      outputBinding:
        glob: '*.csv'
        outputEval: $(inheritMetadata(self, inputs.in_bqsr_reports))
      sbg:fileTypes: CSV

    baseCommand:
    - /opt/gatk-4.1.0.0/gatk --java-options
    arguments:
    - prefix: ''
      position: 1
      valueFrom: |-
        ${
            if (inputs.mem_per_job) {
                return '\"-Xmx'.concat(inputs.mem_per_job, 'M') + '\"';
            } else {
                return '\"-Xmx2048M\"';
            }
        }
      shellQuote: false
    - prefix: ''
      position: 2
      valueFrom: GatherBQSRReports
      shellQuote: false
    - prefix: --output
      position: 3
      valueFrom: |-
        ${
            //sort list of input files by nameroot
            function sortNameroot(x, y) {
                if (x.nameroot < y.nameroot) {
                    return -1;
                }
                if (x.nameroot > y.nameroot) {
                    return 1;
                }
                return 0;
            }
                
            var output_prefix;
            var output_ext;
            var in_num = [].concat(inputs.in_bqsr_reports).length;
            var output_ext = ".csv";
            var in_bqsr_reports = [].concat(inputs.in_bqsr_reports);
            
            //if input_prefix is provided by the user
            if (inputs.prefix) {
                output_prefix = inputs.prefix;
                if (in_num > 1) {
                    output_prefix = output_prefix + '.' + in_num;
                }
            }
            else {
                //if there is only one input file
                if(in_num == 1){
                    // check if the sample_id metadata value is defined for the input file
                    if(in_bqsr_reports[0].metadata && in_bqsr_reports[0].metadata.sample_id) {
                        output_prefix = in_bqsr_reports[0].metadata.sample_id;
                    // if sample_id is not defined
                    } else {
                        output_prefix = in_bqsr_reports[0].path.split('/').pop().split('.')[0];
                    }
                }
                //if there are more than 1 input files
                //sort list of input file objects alphabetically by file name 
                //take the first element from that list, and generate output file name as if that file is the only file on the input. 
                else if(in_num > 1) {
                    //sort list of input files by nameroot
                    in_bqsr_reports.sort(sortNameroot);
                    //take the first alphabetically sorted file
                    var first_file = in_bqsr_reports[0];
                    //check if the sample_id metadata value is defined for the input file
                    if(first_file.metadata && first_file.metadata.sample_id) {
                        output_prefix = first_file.metadata.sample_id + '.' + in_num;
                    // if sample_id is not defined
                    } else {
                        output_prefix = first_file.path.split('/').pop().split('.')[0] + '.' + in_num;
                    }
                }
            }
            var output_full = output_prefix + ".recal_data" + output_ext;
            return output_full;
        }
      shellQuote: false
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-gatherbqsrreports-4-1-0-0/14
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    - CWL1.0
    sbg:content_hash: a0739e0aa57b81afb0485d881aae41db8b23cce8d2153fc5715a7794c934f0edb
    sbg:contributors:
    - nens
    - uros_sipetic
    - marijeta_slavkovic
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1554810073
    sbg:id: h-0df13143/h-68156126/h-70cfae2d/0
    sbg:image_url:
    sbg:latestRevision: 14
    sbg:license: BSD 3-Clause License
    sbg:links:
    - id: https://www.broadinstitute.org/gatk/index.php
      label: Homepage
    - id: https://github.com/broadinstitute/gatk
      label: Source Code
    - id: |-
        https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
      label: Download
    - id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
      label: Publication
    - id: https://gatk.broadinstitute.org/hc/en-us/articles/360036359192-GatherBQSRReports
      label: Documentation
    sbg:modifiedBy: marijeta_slavkovic
    sbg:modifiedOn: 1603192324
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 14
    sbg:revisionNotes: description edited (usage example, memory in description etc)
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1554810073
      sbg:revision: 0
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/8
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1554894740
      sbg:revision: 1
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/11
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557487015
      sbg:revision: 2
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/13
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557734524
      sbg:revision: 3
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/17
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557744219
      sbg:revision: 4
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/22
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000599
      sbg:revision: 5
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/23
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351550
      sbg:revision: 6
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/24
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558451160
      sbg:revision: 7
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/25
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1593698671
      sbg:revision: 8
      sbg:revisionNotes: New wrapper
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1593699134
      sbg:revision: 9
      sbg:revisionNotes: Description review suggestions added
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1593780288
      sbg:revision: 10
      sbg:revisionNotes: performance benchmarking cost edited
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1594045532
      sbg:revision: 11
      sbg:revisionNotes: naming description - added one sentence
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1594045569
      sbg:revision: 12
      sbg:revisionNotes: ''
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1598131313
      sbg:revision: 13
      sbg:revisionNotes: added [].concat to arrays
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1603192324
      sbg:revision: 14
      sbg:revisionNotes: description edited (usage example, memory in description
        etc)
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: out_gathered_bqsr_report
  sbg:x: 1494.5830078125
  sbg:y: 330
- id: gatk_applybqsr_4_1_0_0
  label: GATK ApplyBQSR
  in:
  - id: add_output_sam_program_record
    default: 'true'
  - id: bqsr_recal_file
    source: gatk_gatherbqsrreports_4_1_0_0/out_gathered_bqsr_report
  - id: in_alignments
    source:
    - gatk_setnmmdanduqtags_4_1_0_0/out_alignments
  - id: include_intervals_file
    source: sbg_lines_to_interval_list_abr/out_intervals
  - id: in_reference
    source: in_reference
  - id: static_quantized_quals
    default:
    - 10
    - 20
    - 30
  - id: use_original_qualities
    default: true
  scatter:
  - include_intervals_file
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK ApplyBQSR CWL1.0
    doc: |-
      The **GATK ApplyBQSR** tool recalibrates the base quality scores of an input BAM or CRAM file containing reads.

      This tool performs the second pass in a two-stage process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the base qualities of the input reads based on the recalibration table produced by the **GATK BaseRecalibrator** tool. The goal of this procedure is to correct systematic biases that affect the assignment of base quality scores by the sequencer. The first pass consists of calculating the error empirically and finding patterns in how the error varies with basecall features over all bases. The relevant observations are written to the recalibration table. The second pass consists of applying numerical corrections to each individual basecall, based on the patterns identified in the first step (recorded in the recalibration table), and writing out the recalibrated data to a new BAM or CRAM file [1].

      *A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*

      ###Common Use Cases

      * The **GATK ApplyBQSR** tool requires a BAM or CRAM file on its **Input alignments** (`--input`) input and the covariates table (= recalibration file) generated by the **BaseRecalibrator** tool on its **BQSR recal file** input (`--bqsr-recal-file`). If the input alignments file is in the CRAM format, the reference sequence is required on the **Reference** (`--reference`) input of the tool. The tool generates a new alignments file which contains recalibrated read data on its **Output recalibrated alignments** output.

      * Usage example

      ```
       gatk --java-options "-Xmx2048M" ApplyBQSR \
         --reference reference.fasta \
         --input input.bam \
         --bqsr-recal-file recalibration.table \
         --output output.bam

      ```

      * Original qualities can be retained in the output file under the "OQ" tag if desired. See the **Emit original quals** (`--emit-original-quals`) argument for details [1].

      ###Changes Introduced by Seven Bridges

      * The output file will be prefixed using the **Output name prefix** parameter. If this value is not set, the output name will be generated based on the **Sample ID** metadata value from the input alignments file. If the **Sample ID** value is not set, the name will be inherited from the input alignments file name. In case there are multiple files on the **Input alignments** input, the files will be sorted by name and output file name will be generated based on the first file in the sorted file list, following the rules defined in the previous case. Moreover,  **recalibrated** will be added before the extension of the output file name.

      * The user has a possibility to specify the output file format using the **Output file format** argument. Otherwise, the output file format will be the same as the format of the input file.

      * **Include intervals** (`--intervals`) option is divided into **Include intervals string** and **Include intervals file** options.

      * **Exclude intervals** (`--exclude-intervals`) option is divided into **Exclude intervals string** and **Exclude intervals file** options.

      * The following GATK parameters were excluded from the tool wrapper: `--add-output-vcf-command-line`, `--arguments_file`, `--cloud-index-prefetch-buffer`, `--cloud-prefetch-buffer`, `--create-output-bam-md5`, `--create-output-variant-index`, `--create-output-variant-md5`, `--gatk-config-file`, `--gcs-max-retries`, `--gcs-project-for-requester-pays`, `--help`, `--lenient`, `--QUIET`, `--sites-only-vcf-output`, `--showHidden`, `--tmp-dir`, `--use-jdk-deflater`, `--use-jdk-inflater`, `--verbosity`, `--version`

      ###Common Issues and Important Notes

      *  **Memory per job** (`mem_per_job`) input allows a user to set the desired memory requirement when running a tool or adding it to a workflow. This input should be defined in MB. It is propagated to the Memory requirements part and “-Xmx” parameter of the tool. The default value is 2048MB.
      * **Memory overhead per job** (`mem_overhead_per_job`) input allows a user to set the desired overhead memory when running a tool or adding it to a workflow. This input should be defined in MB. This amount will be added to the Memory per job in the Memory requirements section but it will not be added to the “-Xmx” parameter. The default value is 100MB. 
      * Note: GATK tools that take in mapped read data expect a BAM file as the primary format [2]. More on GATK requirements for mapped sequence data formats can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats).
      * Note: **Input alignments** should have corresponding index files in the same folder. 
      * Note: **Reference** FASTA file should have corresponding .fai (FASTA index) and .dict (FASTA dictionary) files in the same folder. 
      * Note: This tool replaces the use of PrintReads for the application of base quality score recalibration as practiced in earlier versions of GATK (2.x and 3.x) [1].
      * Note: You should only run **ApplyBQSR** with the covariates table created from the input BAM or CRAM file [1].
      * Note: This **Read Filter** (`--read-filter`) is automatically applied to the data by the Engine before processing by **ApplyBQSR** [1]: **WellformedReadFilter**
      * Note: If the **Read filter** (`--read-filter`) option is set to "LibraryReadFilter", the **Library** (`--library`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "PlatformReadFilter", the **Platform filter name** (`--platform-filter-name`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to"PlatformUnitReadFilter", the **Black listed lanes** (`--black-listed-lanes`) option must be set to some value. 
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadGroupBlackListReadFilter", the **Read group black list** (`--read-group-black-list`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadGroupReadFilter", the **Keep read group** (`--keep-read-group`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadLengthReadFilter", the **Max read length** (`--max-read-length`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadNameReadFilter", the **Read name** (`--read-name`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "ReadStrandFilter", the **Keep reverse strand only** (`--keep-reverse-strand-only`) option must be set to some value.
      * Note: If the **Read filter** (`--read-filter`) option is set to "SampleReadFilter", the **Sample** (`--sample`) option must be set to some value.
      * Note: The following options are valid only if an appropriate **Read filter** (`--read-filter`) is specified: **Ambig filter bases** (`--ambig-filter-bases`), **Ambig filter frac** (`--ambig-filter-frac`), **Max fragment length** (`--max-fragment-length`), **Maximum mapping quality** (`--maximum-mapping-quality`), **Minimum mapping quality** (`--minimum-mapping-quality`),  **Do not require soft clips** (`--dont-require-soft-clips-both-ends`), **Filter too short** (`--filter-too-short`), **Min read length** (`--min-read-length`). See the description of each parameter for information on the associated **Read filter**.
      * Note: The wrapper has not been tested for the SAM file type on the **Input alignments** input port.

      ###Performance Benchmarking

      Below is a table describing runtimes and task costs of **GATK ApplyBQSR** for a couple of different samples, executed on the AWS cloud instances:

      | Experiment type |  Input size | Duration |  Cost (on-demand) | Instance (AWS) | 
      |:--------------:|:------------:|:--------:|:-------:|:---------:|
      |     RNA-Seq     |  2.2 GB |   8min   | ~0.07$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     |  6.6 GB |   23min   | ~0.21$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 11 GB |  37min  | ~0.33$ | c4.2xlarge (8 CPUs) | 
      |     RNA-Seq     | 22 GB |  1h 16min  | ~0.68$ | c4.2xlarge (8 CPUs) |

      *Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*

      ###References

      [1] [GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360036725911-ApplyBQSR)

      [2] [GATK Mapped sequence data formats](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: "${\n    return inputs.cpu_per_job ? inputs.cpu_per_job : 1;\n}"
      ramMin: |-
        ${
          var memory = 2048;
          
          if(inputs.mem_per_job) {
          	 memory = inputs.mem_per_job;
          }
          if(inputs.mem_overhead_per_job) {
        	memory += inputs.mem_overhead_per_job;
          }
          else {
             memory += 100;
          }
          return memory;
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/marijeta_slavkovic/gatk-4-1-0-0:0
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib:
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file)) {
                file['metadata'] = {}
            }
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
            return file
        };
        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!o2) {
                return o1;
            };
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
                for (var key in commonMetadata) {
                    if (!(key in example)) {
                        delete commonMetadata[key]
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
                if (o1.secondaryFiles) {
                    o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)
                }
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                    if (o1[i].secondaryFiles) {
                        o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)
                    }
                }
            }
            return o1;
        };

    inputs:
    - id: add_output_sam_program_record
      label: Add output SAM program record
      doc: If true, adds a PG tag to created SAM/BAM/CRAM files.
      type:
      - 'null'
      - name: add_output_sam_program_record
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --add-output-sam-program-record
        position: 4
        shellQuote: false
      sbg:altPrefix: -add-output-sam-program-record
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: ambig_filter_bases
      label: Ambig filter bases
      doc: |-
        Valid only if "AmbiguousBaseReadFilter" is specified:
        Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction. Cannot be used in conjuction with argument(s) ambig-filter-frac.
      type: int?
      inputBinding:
        prefix: --ambig-filter-bases
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: 'null'
    - id: ambig_filter_frac
      label: Ambig filter frac
      doc: |-
        Valid only if "AmbiguousBaseReadFilter" is specified:
        Threshold fraction of ambiguous bases. Cannot be used in conjuction with argument(s) ambig-filter-bases.
      type: float?
      inputBinding:
        prefix: --ambig-filter-frac
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '0.05'
    - id: black_listed_lanes
      label: Black listed lanes
      doc: |-
        Valid only if "PlatformUnitReadFilter" is specified:
        Platform unit (PU) to filter out. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.black_listed_lanes)
              {
                  var bl_lanes = [].concat(inputs.black_listed_lanes);
                  var cmd = [];
                  for (var i = 0; i < bl_lanes.length; i++) 
                  {
                      cmd.push('--black-listed-lanes', bl_lanes[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: bqsr_recal_file
      label: BQSR recal file
      doc: Input recalibration table for BQSR.
      type: File
      inputBinding:
        prefix: --bqsr-recal-file
        position: 4
        shellQuote: false
      sbg:altPrefix: -bqsr
      sbg:category: Required Arguments
      sbg:fileTypes: CSV
    - id: create_output_bam_index
      label: Create output BAM/CRAM index
      doc: If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM
        file.
      type:
      - 'null'
      - name: create_output_bam_index
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --create-output-bam-index
        position: 4
        shellQuote: false
      sbg:altPrefix: -OBI
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
    - id: disable_sequence_dictionary_validation
      label: Disable sequence dictionary validation
      doc: |-
        If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!
      type: boolean?
      inputBinding:
        prefix: --disable-sequence-dictionary-validation
        position: 4
        shellQuote: false
      sbg:altPrefix: -disable-sequence-dictionary-validation
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: disable_tool_default_read_filters
      label: Disable tool default read filters
      doc: |-
        Disable all tool default read filters (warning: many tools will not function correctly without their default read filters on).
      type: boolean?
      inputBinding:
        prefix: --disable-tool-default-read-filters
        position: 4
        shellQuote: false
      sbg:altPrefix: -disable-tool-default-read-filters
      sbg:category: Advanced Arguments
      sbg:toolDefaultValue: 'false'
    - id: dont_require_soft_clips_both_ends
      label: Dont require soft clips both ends
      doc: |-
        Valid only if "OverclippedReadFilter" is specified:
        Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block.
      type: boolean?
      inputBinding:
        prefix: --dont-require-soft-clips-both-ends
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: 'false'
    - id: emit_original_quals
      label: Emit original quals
      doc: Emit original base qualities under the OQ tag.
      type: boolean?
      inputBinding:
        prefix: --emit-original-quals
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: exclude_intervals_file
      label: Exclude intervals file
      doc: One or more genomic intervals to exclude from processing.
      type: File?
      inputBinding:
        prefix: --exclude-intervals
        position: 4
        shellQuote: false
      sbg:altPrefix: -XL
      sbg:category: Optional Arguments
      sbg:fileTypes: BED, LIST, INTERVAL_LIST
      sbg:toolDefaultValue: 'null'
    - id: exclude_intervals_string
      label: Exclude intervals string
      doc: |-
        One or more genomic intervals to exclude from processing. This argument may be specified 0 or more times.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |
          ${
              if (inputs.exclude_intervals_string)
              {
                  var exclude_string = [].concat(inputs.exclude_intervals_string);
                  var cmd = [];
                  for (var i = 0; i < exclude_string.length; i++) 
                  {
                      cmd.push('--exclude-intervals', exclude_string[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -XL
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: filter_too_short
      label: Filter too short
      doc: |-
        Valid only if "OverclippedReadFilter" is specified:
        Minimum number of aligned bases.
      type: int?
      inputBinding:
        prefix: --filter-too-short
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '30'
    - id: global_qscore_prior
      label: Global Qscore prior
      doc: Global Qscore Bayesian prior to use for BQSR.
      type: float?
      inputBinding:
        prefix: --global-qscore-prior
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '-1.0'
    - id: in_alignments
      label: Input alignments
      doc: |-
        BAM/SAM/CRAM file containing reads. This argument must be specified at least once.
      type: File[]
      secondaryFiles:
      - |-
        ${
            var in_alignments = self;
            if (in_alignments.nameext == '.bam' || in_alignments.nameext == '.BAM') {
                return [in_alignments.basename + ".bai", in_alignments.nameroot + ".bai"];
            }
            else if (in_alignments.nameext == ".cram" || in_alignments.nameext == '.CRAM') {
                return [in_alignments.basename + ".crai", in_alignments.nameroot + ".crai", in_alignments.basename + ".bai"];     
            }
            return '';
        }
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.in_alignments) {
                  var alignments = [].concat(inputs.in_alignments);
                  var cmd = [];
                  for (var i=0; i<alignments.length; i++) {
                      cmd.push('--input', alignments[i].path);
                  }
                  return cmd.join(' ');
              } 
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, CRAM
    - id: interval_exclusion_padding
      label: Interval exclusion padding
      doc: Amount of padding (in bp) to add to each interval you are excluding.
      type: int?
      inputBinding:
        prefix: --interval-exclusion-padding
        position: 4
        shellQuote: false
      sbg:altPrefix: -ixp
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: interval_merging_rule
      label: Interval merging rule
      doc: Interval merging rule for abutting intervals.
      type:
      - 'null'
      - name: interval_merging_rule
        type: enum
        symbols:
        - ALL
        - OVERLAPPING_ONLY
      inputBinding:
        prefix: --interval-merging-rule
        position: 4
        shellQuote: false
      sbg:altPrefix: -imr
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: ALL
    - id: interval_padding
      label: Interval padding
      doc: Amount of padding (in bp) to add to each interval you are including.
      type: int?
      inputBinding:
        prefix: --interval-padding
        position: 4
        shellQuote: false
      sbg:altPrefix: -ip
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: interval_set_rule
      label: Interval set rule
      doc: Set merging approach to use for combining interval inputs.
      type:
      - 'null'
      - name: interval_set_rule
        type: enum
        symbols:
        - UNION
        - INTERSECTION
      inputBinding:
        prefix: --interval-set-rule
        position: 4
        shellQuote: false
      sbg:altPrefix: -isr
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: UNION
    - id: include_intervals_file
      label: Include intervals file
      doc: One or more genomic intervals over which to operate.
      type: File?
      inputBinding:
        prefix: --intervals
        position: 4
        shellQuote: false
      sbg:altPrefix: -L
      sbg:category: Optional Arguments
      sbg:fileTypes: BED, LIST, INTERVAL_LIST
      sbg:toolDefaultValue: 'null'
    - id: include_intervals_string
      label: Include intervals string
      doc: |-
        One or more genomic intervals over which to operate. This argument may be specified 0 or more times.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.include_intervals_string)
              {
                  var include_string = [].concat(inputs.include_intervals_string);
                  var cmd = [];
                  for (var i = 0; i < include_string.length; i++) 
                  {
                      cmd.push('--intervals', include_string[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -L
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: keep_read_group
      label: Keep read group
      doc: |-
        Valid only if "ReadGroupReadFilter" is specified:
        The name of the read group to keep. Required.
      type: string?
      inputBinding:
        prefix: --keep-read-group
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: keep_reverse_strand_only
      label: Keep reverse strand only
      doc: |-
        Valid only if "ReadStrandFilter" is specified:
        Keep only reads on the reverse strand. Required.
      type:
      - 'null'
      - name: keep_reverse_strand_only
        type: enum
        symbols:
        - 'true'
        - 'false'
      inputBinding:
        prefix: --keep-reverse-strand-only
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: library
      label: Library
      doc: |-
        Valid only if "LibraryReadFilter" is specified:
        Name of the library to keep. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.library)
              {
                  var lib = [].concat(inputs.library);
                  var cmd = [];
                  for (var i = 0; i < lib.length; i++) 
                  {
                      cmd.push('--library', lib[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -library
      sbg:category: Conditional Arguments for readFilter
    - id: max_fragment_length
      label: Max fragment length
      doc: |-
        Valid only if "FragmentLengthReadFilter" is specified:
        Maximum length of fragment (insert size).
      type: int?
      inputBinding:
        prefix: --max-fragment-length
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '1000000'
    - id: max_read_length
      label: Max read length
      doc: |-
        Valid only if "ReadLengthReadFilter" is specified:
        Keep only reads with length at most equal to the specified value. Required.
      type: int?
      inputBinding:
        prefix: --max-read-length
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: maximum_mapping_quality
      label: Maximum mapping quality
      doc: |-
        Valid only if "MappingQualityReadFilter" is specified:
        Maximum mapping quality to keep (inclusive).
      type: int?
      inputBinding:
        prefix: --maximum-mapping-quality
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: 'null'
    - id: mem_overhead_per_job
      label: Memory overhead per job
      doc: |-
        It allows a user to set the desired overhead memory (in MB) when running a tool or adding it to a workflow.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '100'
    - id: mem_per_job
      label: Memory per job
      doc: |-
        It allows a user to set the desired memory requirement (in MB) when running a tool or adding it to a workflow.in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '2048'
    - id: min_read_length
      label: Min read length
      doc: |-
        Valid only if "ReadLengthReadFilter" is specified:
        Keep only reads with length at least equal to the specified value.
      type: int?
      inputBinding:
        prefix: --min-read-length
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '1'
    - id: minimum_mapping_quality
      label: Minimum mapping quality
      doc: |-
        Valid only if "MappingQualityReadFilter" is specified:
        Minimum mapping quality to keep (inclusive).
      type: int?
      inputBinding:
        prefix: --minimum-mapping-quality
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
      sbg:toolDefaultValue: '10'
    - id: platform_filter_name
      label: Platform filter name
      doc: |-
        Valid only if "PlatformReadFilter" is specified:
        Platform attribute (PL) to match. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.platform_filter_name)
              {
                  var pfn = [].concat(inputs.platform_filter_name);
                  var cmd = [];
                  for (var i = 0; i < pfn.length; i++) 
                  {
                      cmd.push('--platform-filter-name', pfn[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: preserve_qscores_less_than
      label: Preserve qscores less than
      doc: Don't recalibrate bases with quality scores less than this threshold.
      type: int?
      inputBinding:
        prefix: --preserve-qscores-less-than
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '6'
    - id: quantize_quals
      label: Quantize quals
      doc: |-
        Quantize quality scores to a given number of levels. A value of 0 here means "do not quantize". Any value greater than zero will be used to recalculate the quantization using that many levels. Negative values mean that we should quantize using the recalibration report's quantization level. Cannot be used in conjuction with argument(s) static-quantized-quals, round-down-quantized.
      type: int?
      inputBinding:
        prefix: --quantize-quals
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
    - id: read_filter
      label: Read filter
      doc: |-
        Read filters to be applied before analysis. This argument may be specified 0 or more times.
      type:
      - 'null'
      - type: array
        items:
          name: read_filter
          type: enum
          symbols:
          - AlignmentAgreesWithHeaderReadFilter
          - AllowAllReadsReadFilter
          - AmbiguousBaseReadFilter
          - CigarContainsNoNOperator
          - FirstOfPairReadFilter
          - FragmentLengthReadFilter
          - GoodCigarReadFilter
          - HasReadGroupReadFilter
          - LibraryReadFilter
          - MappedReadFilter
          - MappingQualityAvailableReadFilter
          - MappingQualityNotZeroReadFilter
          - MappingQualityReadFilter
          - MatchingBasesAndQualsReadFilter
          - MateDifferentStrandReadFilter
          - MateOnSameContigOrNoMappedMateReadFilter
          - MetricsReadFilter
          - NonChimericOriginalAlignmentReadFilter
          - NonZeroFragmentLengthReadFilter
          - NonZeroReferenceLengthAlignmentReadFilter
          - NotDuplicateReadFilter
          - NotOpticalDuplicateReadFilter
          - NotSecondaryAlignmentReadFilter
          - NotSupplementaryAlignmentReadFilter
          - OverclippedReadFilter
          - PairedReadFilter
          - PassesVendorQualityCheckReadFilter
          - PlatformReadFilter
          - PlatformUnitReadFilter
          - PrimaryLineReadFilter
          - ProperlyPairedReadFilter
          - ReadGroupBlackListReadFilter
          - ReadGroupReadFilter
          - ReadLengthEqualsCigarLengthReadFilter
          - ReadLengthReadFilter
          - ReadNameReadFilter
          - ReadStrandFilter
          - SampleReadFilter
          - SecondOfPairReadFilter
          - SeqIsStoredReadFilter
          - ValidAlignmentEndReadFilter
          - ValidAlignmentStartReadFilter
          - WellformedReadFilter
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--read-filter', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:altPrefix: -RF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: read_group_black_list
      label: Read group black list
      doc: |-
        Valid only if "ReadGroupBlackListReadFilter" is specified:
        The name of the read group to filter out. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.read_group_black_list)
              {
                  var rgbl = [].concat(inputs.read_group_black_list);
                  var cmd = [];
                  for (var i = 0; i < rgbl.length; i++) 
                  {
                      cmd.push('--read-group-black-list', rgbl[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: read_name
      label: Read name
      doc: |-
        Valid only if "ReadNameReadFilter" is specified:
        Keep only reads with this read name. Required.
      type: string?
      inputBinding:
        prefix: --read-name
        position: 4
        shellQuote: false
      sbg:category: Conditional Arguments for readFilter
    - id: read_validation_stringency
      label: Read validation stringency
      doc: |-
        Validation stringency for all SAM/BAM/CRAM files read by this program. The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
      type:
      - 'null'
      - name: read_validation_stringency
        type: enum
        symbols:
        - STRICT
        - LENIENT
        - SILENT
      inputBinding:
        prefix: --read-validation-stringency
        position: 4
        shellQuote: false
      sbg:altPrefix: -VS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: SILENT
    - id: in_reference
      label: Reference
      doc: Reference sequence.
      type: File?
      secondaryFiles:
      - .fai
      - ^.dict
      inputBinding:
        prefix: --reference
        position: 4
        shellQuote: false
      sbg:altPrefix: -R
      sbg:category: Optional Arguments
      sbg:fileTypes: FASTA, FA
      sbg:toolDefaultValue: 'null'
    - id: round_down_quantized
      label: Round down quantized
      doc: |-
        Round quals down to nearest quantized qual. Cannot be used in conjuction with argument quantize-quals.
      type: boolean?
      inputBinding:
        prefix: --round-down-quantized
        position: 4
        shellQuote: false
      sbg:category: Advanced Arguments
      sbg:toolDefaultValue: 'false'
    - id: sample
      label: Sample
      doc: |-
        Valid only if "SampleReadFilter" is specified:
        The name of the sample(s) to keep, filtering out all others. This argument must be specified at least once. Required.
      type: string[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.sample)
              {
                  var samp = [].concat(inputs.sample);
                  var cmd = [];
                  for (var i = 0; i < samp.length; i++) 
                  {
                      cmd.push('--sample', samp[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -sample
      sbg:category: Conditional Arguments for readFilter
    - id: sequence_dictionary
      label: Sequence dictionary
      doc: |-
        Use the given sequence dictionary as the master/canonical sequence dictionary. Must be a .dict file.
      type: File?
      inputBinding:
        prefix: --sequence-dictionary
        position: 4
        shellQuote: false
      sbg:altPrefix: -sequence-dictionary
      sbg:category: Optional Arguments
      sbg:fileTypes: DICT
      sbg:toolDefaultValue: '10.0'
    - id: static_quantized_quals
      label: Static quantized quals
      doc: |-
        Use static quantized quality scores to a given number of levels (with -bqsr). Cannot be used in conjuction with argument(s) quantize-quals. This argument may be specified 0 or more times.
      type: int[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.static_quantized_quals)
              {
                  var sqq = [].concat(inputs.static_quantized_quals);
                  var cmd = [];
                  for (var i = 0; i < sqq.length; i++) 
                  {
                      cmd.push('--static-quantized-quals', sqq[i]);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:category: Advanced Arguments
      sbg:toolDefaultValue: 'null'
    - id: use_original_qualities
      label: Use original qualities
      doc: Use the base quality scores from the OQ tag.
      type: boolean?
      inputBinding:
        prefix: --use-original-qualities
        position: 4
        shellQuote: false
      sbg:altPrefix: -OQ
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: cpu_per_job
      label: CPU per job
      doc: CPU per job.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
    - id: prefix
      label: Output name prefix
      doc: Output file name prefix.
      type: string?
      sbg:category: Config Inputs
    - id: output_extension
      label: Output file format
      doc: Output file format
      type:
      - 'null'
      - name: output_extension
        type: enum
        symbols:
        - sam
        - bam
        - cram
      sbg:category: Config Inputs
    - id: disable_bam_index_caching
      label: Disable BAM index caching
      doc: |-
        If true, don't cache BAM indexes, this will reduce memory requirements but may harm performance if many intervals are specified. Caching is automatically disabled if there are no intervals specified.
      type: boolean?
      inputBinding:
        prefix: --disable-bam-index-caching
        position: 4
        shellQuote: false
      sbg:altPrefix: -DBIC
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: read_index
      label: Read index
      doc: |-
        Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.
      type: File[]?
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (inputs.read_index)
              {
                  var r_index = [].concat(inputs.read_index);
                  var cmd = [];
                  for (var i = 0; i < r_index.length; i++) 
                  {
                      cmd.push('--read-index', r_index[i].path);
                  }
                  return cmd.join(' ');
              }
              return '';
          }
        shellQuote: false
      sbg:altPrefix: -read-index
      sbg:category: Optional Arguments
      sbg:fileTypes: BAI, CRAI
    - id: seconds_between_progress_updates
      label: Seconds between progress updates
      doc: Output traversal statistics every time this many seconds elapse.
      type: float?
      inputBinding:
        prefix: --seconds-between-progress-updates
        position: 4
        shellQuote: false
      sbg:altPrefix: -seconds-between-progress-updates
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '10.00'
    - id: disable_read_filter
      label: Disable read filter
      doc: |-
        Read filters to be disabled before analysis. This argument may be specified 0 or more times.
      type:
      - 'null'
      - type: array
        items:
          name: disable_read_filter
          type: enum
          symbols:
          - WellformedReadFilter
      inputBinding:
        position: 4
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('--disable-read-filter', self[i]);
                  }
                  return cmd.join(' ');
              }
              
          }
        shellQuote: false
      sbg:altPrefix: -DF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'

    outputs:
    - id: out_alignments
      label: Output recalibrated alignments
      doc: Output recalibrated BAM/SAM/CRAM file.
      type: File?
      secondaryFiles:
      - |-
        ${

            if (self.nameext == '.bam' || self.nameext == '.BAM')
            {
                return self.nameroot + ".bai";
            }
            else if (self.nameext == ".cram" || self.nameext == '.CRAM')
            {
                return self.basename + ".bai";     
            }
        }
      outputBinding:
        glob: '*am'
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: BAM, SAM, CRAM

    baseCommand:
    - /opt/gatk-4.1.0.0/gatk --java-options
    arguments:
    - prefix: ''
      position: 1
      valueFrom: |-
        ${
            if (inputs.mem_per_job) {
                return '\"-Xmx'.concat(inputs.mem_per_job, 'M') + '\"';
            }  else {
                return '\"-Xmx2048M\"';
            }
        }
      shellQuote: false
    - prefix: ''
      position: 2
      valueFrom: ApplyBQSR
      shellQuote: false
    - prefix: --output
      position: 3
      valueFrom: |-
        ${
            //sort list of input files by nameroot
            function sortNameroot(x, y) {
                if (x.nameroot < y.nameroot) {
                    return -1;
                }
                if (x.nameroot > y.nameroot) {
                    return 1;
                }
                return 0;
            }
                
            var output_prefix;
            var output_ext;
            var in_num = [].concat(inputs.in_alignments).length;
            var in_align = [].concat(inputs.in_alignments);
            
            //if input_prefix is provided by the user
            if (inputs.prefix) {
                output_prefix = inputs.prefix;
                if (in_num > 1) {
                    output_prefix = output_prefix + '.' + in_num;
                }
            }
            else {
                //if there is only one input file
                if(in_num == 1){
                    // check if the sample_id metadata value is defined for the input file
                    if(in_align[0].metadata && in_align[0].metadata.sample_id) {
                        output_prefix = in_align[0].metadata.sample_id;
                    // if sample_id is not defined
                    } else {
                        output_prefix = in_align[0].path.split('/').pop().split('.')[0];
                    }
                }
                //if there are more than 1 input files
                //sort list of input file objects alphabetically by file name 
                //take the first element from that list, and generate output file name as if that file is the only file on the input. 
                else if(in_num > 1) {
                    //sort list of input files by nameroot
                    in_align.sort(sortNameroot);
                    //take the first alphabetically sorted file
                    var first_file = in_align[0];
                    //check if the sample_id metadata value is defined for the input file
                    if(first_file.metadata && first_file.metadata.sample_id) {
                        output_prefix = first_file.metadata.sample_id + '.' + in_num;
                    // if sample_id is not defined
                    } else {
                        output_prefix = first_file.path.split('/').pop().split('.')[0] + '.' + in_num;
                    }
                }
            }
            output_ext = inputs.output_extension ? inputs.output_extension : in_align[0].path.split('.').pop();
            var output_full = output_prefix + '.recalibrated.' + output_ext;
            return output_full;
        }
      shellQuote: false
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-applybqsr-4-1-0-0/20
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    - CWL1.0
    sbg:content_hash: ae6013c26b9a6948fd717a2ab74f3f08e052e1c3494d04be4ea45b62c71ae729d
    sbg:contributors:
    - uros_sipetic
    - marijeta_slavkovic
    - nens
    - veliborka_josipovic
    - nemanja.vucic
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552923344
    sbg:id: h-11ae251d/h-f4ccbabc/h-67180121/0
    sbg:image_url:
    sbg:latestRevision: 20
    sbg:license: BSD 3-Clause License
    sbg:links:
    - id: https://www.broadinstitute.org/gatk/index.php
      label: Homepage
    - id: https://github.com/broadinstitute/gatk
      label: Source Code
    - id: |-
        https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
      label: Download
    - id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
      label: Publication
    - id: https://gatk.broadinstitute.org/hc/en-us/articles/360036725911-ApplyBQSR
      label: Documentation
    sbg:modifiedBy: marijeta_slavkovic
    sbg:modifiedOn: 1604413284
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 20
    sbg:revisionNotes: small description
    sbg:revisionsInfo:
    - sbg:modifiedBy: uros_sipetic
      sbg:modifiedOn: 1552923344
      sbg:revision: 0
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/8
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554493022
      sbg:revision: 1
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/14
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554493059
      sbg:revision: 2
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/15
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554720859
      sbg:revision: 3
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/16
    - sbg:modifiedBy: veliborka_josipovic
      sbg:modifiedOn: 1554999197
      sbg:revision: 4
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/17
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557734544
      sbg:revision: 5
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/18
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000590
      sbg:revision: 6
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/19
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351541
      sbg:revision: 7
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/21
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558451164
      sbg:revision: 8
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/22
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558524331
      sbg:revision: 9
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/23
    - sbg:modifiedBy: nemanja.vucic
      sbg:modifiedOn: 1559744828
      sbg:revision: 10
      sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/24
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1593698892
      sbg:revision: 11
      sbg:revisionNotes: New wrapper
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1593699993
      sbg:revision: 12
      sbg:revisionNotes: Description review suggestions added
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1594048059
      sbg:revision: 13
      sbg:revisionNotes: naming description and benchmarking cost review
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1594725390
      sbg:revision: 14
      sbg:revisionNotes: added CRAM and SAM to suggested types for in_alignments
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1594725524
      sbg:revision: 15
      sbg:revisionNotes: removed SAM as file suggestion
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1597669898
      sbg:revision: 16
      sbg:revisionNotes: changed default mem_per_job to 2048
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1598131508
      sbg:revision: 17
      sbg:revisionNotes: added [].concat to arrays
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1603199958
      sbg:revision: 18
      sbg:revisionNotes: description edited (usage example Xmx, memory in description
        etc)
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1603296434
      sbg:revision: 19
      sbg:revisionNotes: small description
    - sbg:modifiedBy: marijeta_slavkovic
      sbg:modifiedOn: 1604413284
      sbg:revision: 20
      sbg:revisionNotes: small description
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: out_alignments
  sbg:x: 1615.560546875
  sbg:y: 207.82618713378906
- id: gatk_gatherbamfiles_4_1_0_0
  label: GATK GatherBamFiles
  in:
  - id: create_index
    default: true
  - id: in_alignments
    source:
    - gatk_applybqsr_4_1_0_0/out_alignments
  - id: create_md5_file
    default: true
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: GATK GatherBamFiles
    doc: |-
      **GATK GatherBamFiles** concatenates one or more BAM files resulted form scattered paralel anaysis. 


      ### Common Use Cases 

      * **GATK GatherBamFiles**  tool performs a rapid "gather" or concatenation on BAM files into single BAM file. This is often needed in operations that have been run in parallel across genomics regions by scattering their execution across computing nodes and cores thus resulting in smaller BAM files.
      * Usage example:
      ```

      java -jar picard.jar GatherBamFiles
            --INPUT=input1.bam
            --INPUT=input2.bam
      ```

      ### Common Issues and Important Notes
      * **GATK GatherBamFiles** assumes that the list of BAM files provided as input are in the order that they should be concatenated and simply links the bodies of the BAM files while retaining the header from the first file. 
      *  Operates by copying the gzip blocks directly for speed but also supports the generation of an MD5 in the output file and the indexing of the output BAM file.
      * This tool only support BAM files. It does not support SAM files.

      ###Changes Intorduced by Seven Bridges
      * Generated output BAM file will be prefixed using the **Output prefix** parameter. In case the **Output prefix** is not provided, the output prefix will be the same as the **Sample ID** metadata from the **Input alignments**, if the **Sample ID** metadata exists. Otherwise, the output prefix will be inferred from the **Input alignments** filename. This way, having identical names of the output files between runs is avoided.
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: "${\n    return inputs.cpu_per_job ? inputs.cpu_per_job : 1;\n}"
      ramMin: |-
        ${
            var memory = 4096;
            if (inputs.memory_per_job) 
            {
                memory = inputs.memory_per_job;
            }
            if (inputs.memory_overhead_per_job)
            {
                memory += inputs.memory_overhead_per_job;
            }
            return memory;
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


        var setMetadata = function(file, metadata) {
            if (!('metadata' in file)) {
                file['metadata'] = {}
            }
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

    inputs:
    - id: memory_overhead_per_job
      label: Memory Overhead Per Job
      doc: Memory overhead which will be allocated for one job.
      type: int?
      sbg:category: Execution
    - id: max_records_in_ram
      label: Max records in ram
      doc: |-
        When writing files that need to be sorted, this will specify the number of records stored in ram before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of ram needed.
      type: int?
      inputBinding:
        prefix: --MAX_RECORDS_IN_RAM
        position: 20
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
    - id: memory_per_job
      label: Memory Per Job
      doc: Memory which will be allocated for execution.
      type: int?
      sbg:category: Execution
    - id: create_index
      label: Create index
      doc: Whether to create a bam index when writing a coordinate-sorted bam file.
      type: boolean?
      inputBinding:
        prefix: --CREATE_INDEX
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
    - id: in_reference
      label: Reference sequence
      doc: Reference sequence file.
      type: File?
      inputBinding:
        prefix: --REFERENCE_SEQUENCE
        position: 7
        shellQuote: false
      sbg:altPrefix: -R
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
    - id: output_prefix
      label: Output prefix
      doc: Name of the output bam file to write to.
      type: string?
      sbg:category: Optional Arguments
    - id: in_alignments
      label: Input alignments
      doc: |-
        Two or more bam files or text files containing lists of bam files (one per line). This argument must be specified at least once.
      type: File[]
      inputBinding:
        position: 3
        valueFrom: |-
          ${
             if (self)
             {
                 var cmd = [];
                 for (var i = 0; i < self.length; i++)
                 {
                     cmd.push('--INPUT', self[i].path);
                 }
                 return cmd.join(' ');
             }

          }
        shellQuote: false
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM
    - id: compression_level
      label: Compression level
      doc: Compression level for all compressed files created (e.g. Bam and vcf).
      type: int?
      inputBinding:
        prefix: --COMPRESSION_LEVEL
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
    - id: validation_stringency
      label: Validation stringency
      doc: |-
        Validation stringency for all sam files read by this program. Setting stringency to silent can improve performance when processing a bam file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
      type:
      - 'null'
      - name: validation_stringency
        type: enum
        symbols:
        - STRICT
        - LENIENT
        - SILENT
      inputBinding:
        prefix: --VALIDATION_STRINGENCY
        position: 4
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
    - id: create_md5_file
      label: Create MD5 file
      doc: Whether to create an MD5 digest for any BAM or FASTQ files created.
      type: boolean?
      inputBinding:
        prefix: --CREATE_MD5_FILE
        position: 5
        shellQuote: false
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'FALSE'
    - id: cpu_per_job
      label: CPU per job
      doc: |-
        This input allows a user to set the desired CPU requirement when running a tool or adding it to a workflow.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'

    outputs:
    - id: out_alignments
      label: Output BAM file
      doc: Output BAM file obtained by merging input BAM files.
      type: File?
      secondaryFiles:
      - |-
        ${
            if (inputs.create_index)
            {
                return [self.basename + ".bai", self.nameroot + ".bai"];
            }
            else {
                return ''; 
            }
        }
      outputBinding:
        glob: '*.bam'
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: BAM
    - id: out_md5
      label: MD5 file
      doc: MD5 ouput BAM file.
      type: File?
      outputBinding:
        glob: '*.md5'
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: MD5

    baseCommand: []
    arguments:
    - position: 0
      valueFrom: /opt/gatk --java-options
      shellQuote: false
    - position: 2
      valueFrom: |-
        ${
            if (inputs.memory_per_job) {
                return '\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\"';
            }
            return '\"-Xmx2048M\"';
        }
      shellQuote: false
    - position: 4
      valueFrom: |-
        ${
            var tmp = [].concat(inputs.in_alignments);
                
            if (inputs.output_prefix) {
                return '-O ' +  inputs.output_prefix + ".bam";
                
            }else if (tmp[0].metadata && tmp[0].metadata.sample_id) {
                
                return '-O ' +  tmp[0].metadata.sample_id + ".bam";
            } else {
                 
                return '-O ' +  tmp[0].path.split('/').pop().split(".")[0] + ".bam";
            }
            
            
        }
      shellQuote: false
    - position: 3
      valueFrom: GatherBamFiles
      shellQuote: false
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-gatherbamfiles-4-1-0-0/9
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    sbg:content_hash: adc3fdd806bf7e70cfd29e650f70e8bdc6477baa1d0dc7ef7792f2f8806bcd064
    sbg:contributors:
    - nens
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/23
    sbg:createdBy: nens
    sbg:createdOn: 1554894822
    sbg:id: h-b7eb95de/h-5f9c84aa/h-f6ef8313/0
    sbg:image_url:
    sbg:latestRevision: 9
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - id: https://software.broadinstitute.org/gatk/
      label: Homepage
    - id: |-
        https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_GatherBamFiles.php
      label: Documentation
    - id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
      label: Publications
    - id: https://github.com/broadinstitute/gatk/
      label: Source
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1558531990
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 9
    sbg:revisionNotes: |-
      Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/23
    sbg:revisionsInfo:
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1554894822
      sbg:revision: 0
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/11
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557734548
      sbg:revision: 1
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/14
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1557914509
      sbg:revision: 2
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/16
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558000604
      sbg:revision: 3
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/17
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558351555
      sbg:revision: 4
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/18
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558451620
      sbg:revision: 5
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/19
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558525775
      sbg:revision: 6
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/20
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558526183
      sbg:revision: 7
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/21
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558528334
      sbg:revision: 8
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/22
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1558531990
      sbg:revision: 9
      sbg:revisionNotes: |-
        Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/23
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  out:
  - id: out_alignments
  - id: out_md5
  sbg:x: 1867.5662841796875
  sbg:y: 208.6806640625
- id: samtools_view_1_9_cwl1_0
  label: Samtools View
  in:
  - id: output_format
    default: BAM
  - id: fast_bam_compression
    default: true
  - id: include_header
    default: false
  - id: in_alignments
    source: bwa_mem_bundle_0_7_15/aligned_reads
  scatter:
  - in_alignments
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: Samtools View CWL1.0
    doc: |-
      **SAMtools View** tool prints all alignments from a SAM, BAM, or CRAM file to an output file in SAM format (headerless). You may specify one or more space-separated region specifications to restrict output to only those alignments which overlap the specified region(s). Use of region specifications requires a coordinate-sorted and indexed input file (in BAM or CRAM format) [1].

      *A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*

      ####Regions

      Regions can be specified as: RNAME[:STARTPOS[-ENDPOS]] and all position coordinates are 1-based. 

      **Important note:** when multiple regions are given, some alignments may be output multiple times if they overlap more than one of the specified regions.

      Examples of region specifications:

      - **chr1**  - Output all alignments mapped to the reference sequence named `chr1' (i.e. @SQ SN:chr1).

      - **chr2:1000000** - The region on chr2 beginning at base position 1,000,000 and ending at the end of the chromosome.

      - **chr3:1000-2000** - The 1001bp region on chr3 beginning at base position 1,000 and ending at base position 2,000 (including both end positions).

      - **'\*'** - Output the unmapped reads at the end of the file. (This does not include any unmapped reads placed on a reference sequence alongside their mapped mates.)

      - **.** - Output all alignments. (Mostly unnecessary as not specifying a region at all has the same effect.) [1]

      ###Common Use Cases

      This tool can be used for: 

      - Filtering BAM/SAM/CRAM files - options set by the following parameters and input files: **Include reads with all of these flags** (`-f`), **Exclude reads with any of these flags** (`-F`), **Exclude reads with all of these flags** (`-G`), **Read group** (`-r`), **Minimum mapping quality** (`-q`), **Only include alignments in library** (`-l`), **Minimum number of CIGAR bases consuming query sequence** (`-m`), **Subsample fraction** (`-s`), **Read group list** (`-R`), **BED region file** (`-L`)
      - Format conversion between SAM/BAM/CRAM formats - set by the following parameters: **Output format** (`--output-fmt/-O`), **Fast bam compression** (`-1`), **Output uncompressed BAM** (`-u`)
      - Modification of the data which is contained in each alignment - set by the following parameters: **Collapse the backward CIGAR operation** (`-B`), **Read tags to strip** (`-x`)
      - Counting number of alignments in SAM/BAM/CRAM file - set by parameter **Output only count of matching records** (`-c`)

      ###Changes Introduced by Seven Bridges

      - Parameters **Output BAM** (`-b`) and **Output CRAM** (`-C`) were excluded from the wrapper since they are redundant with parameter **Output format** (`--output-fmt/-O`).
      - Parameter **Input format** (`-S`) was excluded from wrapper since it is ignored by the tool (input format is auto-detected).
      - Input file **Index file** was added to the wrapper to enable operations that require an index file for BAM/CRAM files.
      - Parameter **Number of threads** (`--threads/-@`) specifies the total number of threads instead of additional threads. Command line argument (`--threads/-@`) will be reduced by 1 to set the number of additional threads.

      ###Common Issues and Important Notes

      - When multiple regions are given, some alignments may be output multiple times if they overlap more than one of the specified regions [1].
      - Use of region specifications requires a coordinate-sorted and indexed input file (in BAM or CRAM format) [1].
      - Option **Output uncompressed BAM** (`-u`) saves time spent on compression/decompression and is thus preferred when the output is piped to another SAMtools command [1].

      ###Performance Benchmarking

      Multithreading can be enabled by setting parameter **Number of threads** (`--threads/-@`). In the following table you can find estimates of **SAMtools View** running time and cost. 

      *Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*  

      | Input type | Input size | # of reads | Read length | Output format | # of threads | Duration | Cost | Instance (AWS)|
      |---------------|--------------|-----------------|---------------|------------------|-------------------|-----------------|-------------|--------|-------------|
      | BAM | 5.26 GB | 71.5M | 76 | BAM | 1 | 13min. | \$0.12 | c4.2xlarge |
      | BAM | 11.86 GB | 161.2M | 101 | BAM | 1 | 33min. | \$0.30 | c4.2xlarge |
      | BAM | 18.36 GB | 179M | 76 | BAM | 1 | 60min. | \$0.54 | c4.2xlarge |
      | BAM | 58.61 GB | 845.6M | 150 | BAM | 1 | 3h 25min. | \$1.84 | c4.2xlarge |
      | BAM | 5.26 GB | 71.5M | 76 | BAM | 8 | 5min. | \$0.04 | c4.2xlarge |
      | BAM | 11.86 GB | 161.2M | 101 | BAM | 8 | 11min. | \$0.10 | c4.2xlarge |
      | BAM | 18.36 GB | 179M | 76 | BAM | 8 | 19min. | \$0.17 | c4.2xlarge |
      | BAM | 58.61 GB | 845.6M | 150 | BAM | 8 | 61min. | \$0.55 | c4.2xlarge |
      | BAM | 5.26 GB | 71.5M | 76 | SAM | 8 | 14min. | \$0.13 | c4.2xlarge |
      | BAM | 11.86 GB | 161.2M | 101 | SAM | 8 | 23min. | \$0.21 | c4.2xlarge |
      | BAM | 18.36 GB | 179M | 76 | SAM | 8 | 35min. | \$0.31 | c4.2xlarge |
      | BAM | 58.61 GB | 845.6M | 150 | SAM | 8 | 2h 29min. | \$1.34 | c4.2xlarge |

      ###References

      [1] [SAMtools documentation](http://www.htslib.org/doc/samtools-1.9.html)
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: |-
        ${
          if (inputs.cpu_per_job) {
              return inputs.cpu_per_job
          }
          else {
          if((inputs.threads)){
            return (inputs.threads)
          }
          else{
            return 1
          }
          }
        }
      ramMin: |-
        ${
          if (inputs.mem_per_job) {
              return inputs.mem_per_job
          }    
          else {
          mem_offset = 1000
          if((inputs.in_reference)){
            mem_offset = mem_offset + 3000
          }
          if((inputs.threads)){
            threads = (inputs.threads)
          }
          else{
            threads = 1
          }
          return mem_offset + threads * 500
          }
        }
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/jrandjelovic/samtools-1-9:1
    - class: InitialWorkDirRequirement
      listing:
      - $(inputs.in_reference)
      - $(inputs.reference_file_list)
      - $(inputs.in_index)
      - $(inputs.in_alignments)
    - class: InlineJavascriptRequirement
      expressionLib:
      - |2-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

    inputs:
    - id: in_index
      label: Index file
      doc: This tool requires index file for some use cases.
      type: File?
      sbg:category: File inputs
      sbg:fileTypes: BAI, CRAI, CSI
    - id: output_format
      label: Output format
      doc: Output file format
      type:
      - 'null'
      - name: output_format
        type: enum
        symbols:
        - SAM
        - BAM
        - CRAM
      inputBinding:
        prefix: --output-fmt
        position: 1
        shellQuote: false
      sbg:altPrefix: -O
      sbg:category: Config inputs
      sbg:toolDefaultValue: SAM
    - id: fast_bam_compression
      label: Fast BAM compression
      doc: Enable fast BAM compression (implies output in bam format).
      type: boolean?
      inputBinding:
        prefix: '-1'
        position: 2
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: 'False'
    - id: uncompressed_bam
      label: Output uncompressed BAM
      doc: |-
        Output uncompressed BAM (implies output BAM format). This option saves time spent on compression/decompression and is thus preferred when the output is piped to another SAMtools command.
      type: boolean?
      inputBinding:
        prefix: -u
        position: 3
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: 'False'
    - id: include_header
      label: Include the header in the output
      doc: Include the header in the output.
      type: boolean?
      inputBinding:
        prefix: -h
        position: 4
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: 'False'
    - id: output_header_only
      label: Output the header only
      doc: Output the header only.
      type: boolean?
      inputBinding:
        prefix: -H
        position: 5
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: 'False'
    - id: collapse_cigar
      label: Collapse the backward CIGAR operation
      doc: Collapse the backward CIGAR operation.
      type: boolean?
      inputBinding:
        prefix: -B
        position: 6
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: 'False'
    - id: filter_include
      label: Include reads with all of these flags
      doc: |-
        Only output alignments with all bits set in this integer present in the FLAG field.
      type: int?
      inputBinding:
        prefix: -f
        position: 7
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: '0'
    - id: filter_exclude_any
      label: Exclude reads with any of these flags
      doc: |-
        Do not output alignments with any bits set in this integer present in the FLAG field.
      type: int?
      inputBinding:
        prefix: -F
        position: 8
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: '0'
    - id: filter_exclude_all
      label: Exclude reads with all of these flags
      doc: |-
        Only exclude reads with all of the bits set in this integer present in the FLAG field.
      type: int?
      inputBinding:
        prefix: -G
        position: 9
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: '0'
    - id: read_group
      label: Read group
      doc: Only output reads in the specified read group.
      type: string?
      inputBinding:
        prefix: -r
        position: 10
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: 'null'
    - id: filter_mapq
      label: Minimum mapping quality
      doc: Skip alignments with MAPQ smaller than this value.
      type: int?
      inputBinding:
        prefix: -q
        position: 11
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: '0'
    - id: filter_library
      label: Only include alignments in library
      doc: Only output alignments in this library.
      type: string?
      inputBinding:
        prefix: -l
        position: 12
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: 'null'
    - id: min_cigar_operations
      label: Minimum number of CIGAR bases consuming query sequence
      doc: |-
        Only output alignments with number of CIGAR bases consuming query sequence  ≥ INT.
      type: int?
      inputBinding:
        prefix: -m
        position: 13
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: '0'
    - id: read_tag_to_strip
      label: Read tags to strip
      doc: Read tag to exclude from output (repeatable).
      type: string[]?
      inputBinding:
        prefix: ''
        position: 14
        valueFrom: |-
          ${
              if (self)
              {
                  var cmd = [];
                  for (var i = 0; i < self.length; i++) 
                  {
                      cmd.push('-x', self[i]);
                      
                  }
                  return cmd.join(' ');
              }
          }
        itemSeparator: ' '
        shellQuote: false
      sbg:category: Config Inputs
    - id: count_alignments
      label: Output only count of matching records
      doc: |-
        Instead of outputing the alignments, only count them and output the total number. All filter options, such as -f, -F, and -q, are taken into account.
      type: boolean?
      inputBinding:
        prefix: -c
        position: 15
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: 'False'
    - id: input_fmt_option
      label: Input file format option
      doc: Specify a single input file format option in the form of OPTION or OPTION=VALUE.
      type: string?
      inputBinding:
        prefix: --input-fmt-option
        position: 16
        shellQuote: false
      sbg:category: Config Inputs
    - id: output_fmt_option
      label: Output file format option
      doc: |-
        Specify a single output file format option in the form of OPTION or OPTION=VALUE.
      type: string?
      inputBinding:
        prefix: --output-fmt-option
        position: 17
        shellQuote: false
      sbg:category: Config Inputs
    - id: subsample_fraction
      label: Subsample fraction
      doc: |-
        Output only a proportion of the input alignments. This subsampling acts in the same way on all of the alignment records in the same template or read pair, so it never keeps a read but not its mate. The integer and fractional parts of the INT.FRAC are used separately: the part after the decimal point sets the fraction of templates/pairs to be kept, while the integer part is used as a seed that influences which subset of reads is kept. When subsampling data that has previously been subsampled, be sure to use a different seed value from those used previously; otherwise more reads will be retained than expected.
      type: float?
      inputBinding:
        prefix: -s
        position: 18
        shellQuote: false
      sbg:category: Config Inputs
    - id: threads
      label: Number of threads
      doc: |-
        Number of threads. SAMtools uses argument --threads/-@ to specify number of additional threads. This parameter sets total number of threads (and CPU cores). Command line argument will be reduced by 1 to set number of additional threads.
      type: int?
      inputBinding:
        prefix: --threads
        position: 19
        valueFrom: |-
          ${
            if((inputs.threads)){
              return (inputs.threads) - 1
            }
            else{
              return
            }
          }
        shellQuote: false
      sbg:altPrefix: -@
      sbg:category: Execution
      sbg:toolDefaultValue: '1'
    - id: omitted_reads_filename
      label: Filename for reads not selected by filters
      doc: |-
        Write alignments that are not selected by the various filter options to this file. When this option is used, all alignments (or all alignments intersecting the regions specified) are written to either the output file or this file, but never both.
      type: string?
      inputBinding:
        prefix: -U
        position: 20
        shellQuote: false
      sbg:category: Config Inputs
    - id: output_filename
      label: Output filename
      doc: Define a filename of the output.
      type: string?
      default: default_output_filename
      inputBinding:
        prefix: -o
        position: 21
        valueFrom: |-
          ${
            if (inputs.output_filename!="default_output_filename"){
              return (inputs.output_filename)
            }
            input_filename = [].concat(inputs.in_alignments)[0].path.split('/').pop()
            input_name_base = input_filename.split('.').slice(0,-1).join('.')
            ext = 'sam'
            if (inputs.count_alignments){
              return input_name_base + '.count.txt'
            }
            if ((inputs.uncompressed_bam) || (inputs.fast_bam_compression)){
              ext = 'bam'
            }
            if (inputs.output_format){
              ext = (inputs.output_format).toLowerCase()
            }
            if (inputs.output_header_only){
              ext = 'header.' + ext
            }
            if (inputs.subsample_fraction){
              ext = 'subsample.' + ext
            }
            if ((inputs.bed_file) || (inputs.read_group) || (inputs.read_group_list) ||
                (inputs.filter_mapq) || (inputs.filter_library) || (inputs.min_cigar_operations) ||
                (inputs.filter_include) || (inputs.filter_exclude_any) || 
                (inputs.filter_exclude_all) || (inputs.regions_array)){
              ext = 'filtered.' + ext
            }
              
            return input_name_base + '.' + ext
          }
        shellQuote: false
      sbg:category: Config Inputs
      sbg:toolDefaultValue: stdout
    - id: bed_file
      label: BED region file
      doc: Only output alignments overlapping the input BED file.
      type: File?
      inputBinding:
        prefix: -L
        position: 22
        shellQuote: false
      sbg:category: File Inputs
      sbg:fileTypes: BED
    - id: read_group_list
      label: Read group list
      doc: Output alignments in read groups listed in this file.
      type: File?
      inputBinding:
        prefix: -R
        position: 23
        shellQuote: false
      sbg:category: File Inputs
      sbg:fileTypes: TXT
    - id: in_reference
      label: Reference file
      doc: |-
        A FASTA format reference file, optionally compressed by bgzip and ideally indexed by SAMtools Faidx. If an index is not present, one will be generated for you. This file is used for compression/decompression of CRAM files. Please provide reference file when using CRAM input/output file.
      type: File?
      inputBinding:
        prefix: --reference
        position: 24
        shellQuote: false
      sbg:altPrefix: -T
      sbg:category: File Inputs
      sbg:fileTypes: FASTA, FA, FASTA.GZ, FA.GZ, GZ
    - id: reference_file_list
      label: List of reference names and lengths
      doc: |-
        A tab-delimited file. Each line must contain the reference name in the first column and the length of the reference in the second column, with one line for each distinct reference. Any additional fields beyond the second column are ignored. This file also defines the order of the reference sequences in sorting. If you run SAMtools Faidx on reference FASTA file (<ref.fa>), the resulting index file <ref.fa>.fai can be used as this file.
      type: File?
      inputBinding:
        prefix: -t
        position: 25
        shellQuote: false
      sbg:category: File Inputs
      sbg:fileTypes: FAI, TSV, TXT
    - id: in_alignments
      label: Input BAM/SAM/CRAM file
      doc: Input BAM/SAM/CRAM file.
      type: File
      inputBinding:
        position: 99
        shellQuote: false
      sbg:category: File Inputs
      sbg:fileTypes: BAM, SAM, CRAM
    - id: regions_array
      label: Regions array
      doc: |-
        With no options or regions specified, prints all alignments in the specified input alignment file (in SAM, BAM, or CRAM format) to output file in specified format. Use of region specifications requires a coordinate-sorted and indexed input file (in BAM or CRAM format). Regions can be specified as: RNAME[:STARTPOS[-ENDPOS]] and all position coordinates are 1-based.  Important note: when multiple regions are given, some alignments may be output multiple times if they overlap more than one of the specified regions. Examples of region specifications:  chr1 - Output all alignments mapped to the reference sequence named `chr1' (i.e. @SQ SN:chr1);  chr2:1000000 - The region on chr2 beginning at base position 1,000,000 and ending at the end of the chromosome;  chr3:1000-2000 - The 1001bp region on chr3 beginning at base position 1,000 and ending at base position 2,000 (including both end positions);  '*' - Output the unmapped reads at the end of the file (this does not include any unmapped reads placed on a reference sequence alongside their mapped mates.);  . - Output all alignments (mostly unnecessary as not specifying a region at all has the same effect).
      type: string[]?
      inputBinding:
        position: 100
        shellQuote: false
      sbg:category: Config Inputs
    - id: multi_region_iterator
      label: Use the multi-region iterator
      doc: |-
        Use the multi-region iterator on the union of the BED file and command-line region arguments.
      type: boolean?
      inputBinding:
        prefix: -M
        position: 22
        shellQuote: false
      sbg:category: Config inputs
      sbg:toolDefaultValue: 'False'
    - id: mem_per_job
      label: Memory per job
      doc: Memory per job in MB.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1500'
    - id: cpu_per_job
      label: CPU per job
      doc: Number of CPUs per job.
      type: int?
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'

    outputs:
    - id: out_alignments
      label: Output BAM, SAM, or CRAM file
      doc: The output file.
      type: File?
      outputBinding:
        glob: |-
          ${
            if ((inputs.output_filename!="default_output_filename")){
              return (inputs.output_filename)
            }
            input_filename = [].concat((inputs.in_alignments))[0].path.split('/').pop()
            input_name_base = input_filename.split('.').slice(0,-1). join('.')
            ext = 'sam'
            if ((inputs.count_alignments)){
              return 
            }
            if ((inputs.uncompressed_bam) || (inputs.fast_bam_compression)){
              ext = 'bam'
            }
            if ((inputs.output_format)){
              ext = (inputs.output_format).toLowerCase()
            }
            if ((inputs.output_header_only)){
              ext = 'header.' + ext
            }
            if ((inputs.subsample_fraction)){
              ext = 'subsample.' + ext
            }
            if ((inputs.bed_file) || (inputs.read_group) || (inputs.read_group_list) ||
                (inputs.filter_mapq) || (inputs.filter_library) || (inputs.min_cigar_operations) ||
                (inputs.filter_include) || (inputs.filter_exclude_any) || 
                (inputs.filter_exclude_all) || (inputs.regions_array)){
              ext = 'filtered.' + ext
            }
              
            return input_name_base + '.' + ext
          }
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: BAM, SAM, CRAM
    - id: reads_not_selected_by_filters
      label: Reads not selected by filters
      doc: File containing reads that are not selected by filters.
      type: File?
      outputBinding:
        glob: |-
          ${
            if ((inputs.omitted_reads_filename)){
              return (inputs.omitted_reads_filename)
            }
          }
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: BAM, SAM, CRAM
    - id: alignement_count
      label: Alignment count
      doc: File containing number of alignments.
      type: File?
      outputBinding:
        glob: |-
          ${
            input_filename = [].concat((inputs.in_alignments))[0].path.split('/').pop()
            input_name_base = input_filename.split('.').slice(0,-1). join('.')
            return input_name_base + '.count.txt'
          }
        outputEval: $(inheritMetadata(self, inputs.in_alignments))
      sbg:fileTypes: TXT

    baseCommand:
    - /opt/samtools-1.9/samtools
    - view
    id: lea_lenhardt_ackovic/samtools-1-9-cwl1-0-demo/samtools-view-1-9-cwl1-0/6
    sbg:appVersion:
    - v1.0
    sbg:categories:
    - Utilities
    - BAM Processing
    - CWL1.0
    sbg:content_hash: ab372090457bac69a1b2bd8deff4ef40ca29052f82dd4850241d8d9e1096eed34
    sbg:contributors:
    - lea_lenhardt_ackovic
    sbg:createdBy: lea_lenhardt_ackovic
    sbg:createdOn: 1572600501
    sbg:id: h-ed1cbd9c/h-0e016058/h-0e3dbed4/0
    sbg:image_url:
    sbg:latestRevision: 6
    sbg:license: MIT License
    sbg:links:
    - id: http://www.htslib.org/
      label: Homepage
    - id: https://github.com/samtools/samtools
      label: Source Code
    - id: https://github.com/samtools/samtools/wiki
      label: Wiki
    - id: https://sourceforge.net/projects/samtools/files/samtools/
      label: Download
    - id: http://www.ncbi.nlm.nih.gov/pubmed/19505943
      label: Publication
    - id: http://www.htslib.org/doc/samtools-1.9.html
      label: Documentation
    sbg:modifiedBy: lea_lenhardt_ackovic
    sbg:modifiedOn: 1578571408
    sbg:project: lea_lenhardt_ackovic/samtools-1-9-cwl1-0-demo
    sbg:projectName: SAMtools 1.9 - CWL1.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 6
    sbg:revisionNotes: Added file requirements for in_index and in_alignments
    sbg:revisionsInfo:
    - sbg:modifiedBy: lea_lenhardt_ackovic
      sbg:modifiedOn: 1572600501
      sbg:revision: 0
      sbg:revisionNotes:
    - sbg:modifiedBy: lea_lenhardt_ackovic
      sbg:modifiedOn: 1572600525
      sbg:revision: 1
      sbg:revisionNotes: Final version
    - sbg:modifiedBy: lea_lenhardt_ackovic
      sbg:modifiedOn: 1575029042
      sbg:revision: 2
      sbg:revisionNotes: Edited description, tag, default values.
    - sbg:modifiedBy: lea_lenhardt_ackovic
      sbg:modifiedOn: 1575042426
      sbg:revision: 3
      sbg:revisionNotes: mem_per_job default value set
    - sbg:modifiedBy: lea_lenhardt_ackovic
      sbg:modifiedOn: 1576241025
      sbg:revision: 4
      sbg:revisionNotes: Description edited - references put before full stop
    - sbg:modifiedBy: lea_lenhardt_ackovic
      sbg:modifiedOn: 1576242427
      sbg:revision: 5
      sbg:revisionNotes: Categories edited
    - sbg:modifiedBy: lea_lenhardt_ackovic
      sbg:modifiedOn: 1578571408
      sbg:revision: 6
      sbg:revisionNotes: Added file requirements for in_index and in_alignments
    sbg:sbgMaintained: false
    sbg:toolAuthor: |-
      Heng Li (Sanger Institute), Bob Handsaker (Broad Institute), Jue Ruan (Beijing Genome Institute), Colin Hercus, Petr Danecek
    sbg:toolkit: samtools
    sbg:toolkitVersion: '1.9'
    sbg:validationErrors: []
  out:
  - id: out_alignments
  - id: reads_not_selected_by_filters
  - id: alignement_count
  sbg:x: -106.09046173095703
  sbg:y: 247.76466369628906
- id: sbg_lines_to_interval_list_abr
  label: SBG Lines to Interval List
  in:
  - id: input_tsv
    source: gatk_createsequencegroupingtsv_4_1_0_0/sequence_grouping_with_unmapped
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: SBG Lines to Interval List
    doc: |-
      This tools is used for splitting GATK sequence grouping file into subgroups.

      ### Common Use Cases

      Each subgroup file contains intervals defined on single line in grouping file. Grouping file is output of GATKs **CreateSequenceGroupingTSV** script which is used in best practice workflows sush as **GATK Best Practice Germline Workflow**.
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 1000
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/uros_sipetic/sci-python:2.7
    - class: InitialWorkDirRequirement
      listing:
      - entryname: lines_to_intervals.py
        writable: false
        entry: |
          import sys
          import hashlib
          import os
          import json

          obj_template = {
              'basename': '',
              'checksum': '',
              'class': 'File',
              'dirname': '',
              'location': '',
              'nameext': 'intervals',
              'nameroot': '',
              'path': '',
              'size': '',
          }

          with open(sys.argv[1], 'r') as f:

              obj_list = []
              sys.stderr.write('Reading file {}\n'.format(sys.argv[1]))
              nameroot = '.'.join(sys.argv[1].split('/')[-1].split('.')[:-1])
              for i, line in enumerate(f):
                  out_file_name = '{}.group.{}.intervals'.format(nameroot, i+1)
                  out_file = open(out_file_name, 'a')
                  for interval in line.split():
                      out_file.write(interval + '\n')
                  out_file.close()
                  sys.stderr.write('Finished writing to file {}\n'.format(out_file_name))

                  obj = dict(obj_template)
                  obj['basename'] = out_file_name
                  obj['checksum'] = 'sha1$' + hashlib.sha1(open(out_file_name, 'r').read()).hexdigest()
                  obj['dirname'] = os.getcwd()
                  obj['location'] = '/'.join([os.getcwd(), out_file_name])
                  obj['nameroot'] = '.'.join(out_file_name.split('.')[:-1])
                  obj['path'] = '/'.join([os.getcwd(), out_file_name])
                  obj['size'] = os.path.getsize('/'.join([os.getcwd(), out_file_name]))

                  obj_list.append(obj)

              out_json = {'out_intervals': obj_list}

              json.dump(out_json, open('cwl.output.json', 'w'), indent=1)
              sys.stderr.write('Job done.\n')
    - class: InlineJavascriptRequirement

    inputs:
    - id: input_tsv
      label: Input group file
      doc: This file is output of GATKs CreateSequenceGroupingTSV script.
      type: File
      inputBinding:
        position: 1
        shellQuote: false
      sbg:category: Required Arguments
      sbg:fileTypes: TSV, TXT

    outputs:
    - id: out_intervals
      label: Intervals
      doc: GATK Intervals files.
      type: File[]
      sbg:fileTypes: INTERVALS, BED

    baseCommand:
    - python
    - lines_to_intervals.py
    id: sevenbridges/sbgtools-cwl1-0-demo/sbg-lines-to-interval-list/3
    sbg:appVersion:
    - v1.0
    sbg:content_hash: a7c4b064a52abdea428818baaba8fdc326902195b3a61fdfdd774c657825c5cc6
    sbg:contributors:
    - nens
    sbg:createdBy: nens
    sbg:createdOn: 1566809066
    sbg:id: h-a73bb3af/h-c5b233b9/h-8f7ccefa/0
    sbg:image_url:
    sbg:latestRevision: 3
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1611663678
    sbg:project: sevenbridges/sbgtools-cwl1-0-demo
    sbg:projectName: SBGTools - CWL1.x - Demo
    sbg:publisher: sbg
    sbg:revision: 3
    sbg:revisionNotes: docker image
    sbg:revisionsInfo:
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1566809066
      sbg:revision: 0
      sbg:revisionNotes:
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1566809311
      sbg:revision: 1
      sbg:revisionNotes: v1 - dev
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1611663319
      sbg:revision: 2
      sbg:revisionNotes: v2 - dev
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1611663678
      sbg:revision: 3
      sbg:revisionNotes: docker image
    sbg:sbgMaintained: false
    sbg:toolAuthor: Stefan Stojanovic
    sbg:toolkit: SBG Tools
    sbg:toolkitVersion: '1.0'
    sbg:validationErrors: []
  out:
  - id: out_intervals
  sbg:x: 981.438232421875
  sbg:y: -67.39484405517578
- id: sbg_lines_to_interval_list_br
  label: SBG Lines to Interval List
  in:
  - id: input_tsv
    source: gatk_createsequencegroupingtsv_4_1_0_0/sequence_grouping
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: SBG Lines to Interval List
    doc: |-
      This tools is used for splitting GATK sequence grouping file into subgroups.

      ### Common Use Cases

      Each subgroup file contains intervals defined on single line in grouping file. Grouping file is output of GATKs **CreateSequenceGroupingTSV** script which is used in best practice workflows sush as **GATK Best Practice Germline Workflow**.
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 1000
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/uros_sipetic/sci-python:2.7
    - class: InitialWorkDirRequirement
      listing:
      - entryname: lines_to_intervals.py
        writable: false
        entry: |
          import sys
          import hashlib
          import os
          import json

          obj_template = {
              'basename': '',
              'checksum': '',
              'class': 'File',
              'dirname': '',
              'location': '',
              'nameext': 'intervals',
              'nameroot': '',
              'path': '',
              'size': '',
          }

          with open(sys.argv[1], 'r') as f:

              obj_list = []
              sys.stderr.write('Reading file {}\n'.format(sys.argv[1]))
              nameroot = '.'.join(sys.argv[1].split('/')[-1].split('.')[:-1])
              for i, line in enumerate(f):
                  out_file_name = '{}.group.{}.intervals'.format(nameroot, i+1)
                  out_file = open(out_file_name, 'a')
                  for interval in line.split():
                      out_file.write(interval + '\n')
                  out_file.close()
                  sys.stderr.write('Finished writing to file {}\n'.format(out_file_name))

                  obj = dict(obj_template)
                  obj['basename'] = out_file_name
                  obj['checksum'] = 'sha1$' + hashlib.sha1(open(out_file_name, 'r').read()).hexdigest()
                  obj['dirname'] = os.getcwd()
                  obj['location'] = '/'.join([os.getcwd(), out_file_name])
                  obj['nameroot'] = '.'.join(out_file_name.split('.')[:-1])
                  obj['path'] = '/'.join([os.getcwd(), out_file_name])
                  obj['size'] = os.path.getsize('/'.join([os.getcwd(), out_file_name]))

                  obj_list.append(obj)

              out_json = {'out_intervals': obj_list}

              json.dump(out_json, open('cwl.output.json', 'w'), indent=1)
              sys.stderr.write('Job done.\n')
    - class: InlineJavascriptRequirement

    inputs:
    - id: input_tsv
      label: Input group file
      doc: This file is output of GATKs CreateSequenceGroupingTSV script.
      type: File
      inputBinding:
        position: 1
        shellQuote: false
      sbg:category: Required Arguments
      sbg:fileTypes: TSV, TXT

    outputs:
    - id: out_intervals
      label: Intervals
      doc: GATK Intervals files.
      type: File[]
      sbg:fileTypes: INTERVALS, BED

    baseCommand:
    - python
    - lines_to_intervals.py
    id: sevenbridges/sbgtools-cwl1-0-demo/sbg-lines-to-interval-list/3
    sbg:appVersion:
    - v1.0
    sbg:content_hash: a7c4b064a52abdea428818baaba8fdc326902195b3a61fdfdd774c657825c5cc6
    sbg:contributors:
    - nens
    sbg:createdBy: nens
    sbg:createdOn: 1566809066
    sbg:id: h-6005657c/h-06d45b7e/h-c5caf58c/0
    sbg:image_url:
    sbg:latestRevision: 3
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1611663678
    sbg:project: sevenbridges/sbgtools-cwl1-0-demo
    sbg:projectName: SBGTools - CWL1.x - Demo
    sbg:publisher: sbg
    sbg:revision: 3
    sbg:revisionNotes: docker image
    sbg:revisionsInfo:
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1566809066
      sbg:revision: 0
      sbg:revisionNotes:
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1566809311
      sbg:revision: 1
      sbg:revisionNotes: v1 - dev
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1611663319
      sbg:revision: 2
      sbg:revisionNotes: v2 - dev
    - sbg:modifiedBy: nens
      sbg:modifiedOn: 1611663678
      sbg:revision: 3
      sbg:revisionNotes: docker image
    sbg:sbgMaintained: false
    sbg:toolAuthor: Stefan Stojanovic
    sbg:toolkit: SBG Tools
    sbg:toolkitVersion: '1.0'
    sbg:validationErrors: []
  out:
  - id: out_intervals
  sbg:x: 979.7381591796875
  sbg:y: 135.31478881835938

hints:
- class: sbg:AWSInstanceType
  value: c5.9xlarge;ebs-gp2;3000
id: |-
  https://cgc-api.sbgenomics.com/v2/apps/admin/sbg-public-data/broad-best-practice-data-pre-processing-workflow-4-1-0-0/26/raw/
sbg:appVersion:
- v1.0
sbg:categories:
- Genomics
- Alignment
- CWL1.0
- GATK
sbg:content_hash: ae9f89d0093c72c279ae8a547a68502dd44cae8f53c8672c8be075476d2962f06
sbg:contributors:
- admin
sbg:createdBy: admin
sbg:createdOn: 1572002743
sbg:expand_workflow: false
sbg:id: |-
  admin/sbg-public-data/broad-best-practice-data-pre-processing-workflow-4-1-0-0/26
sbg:image_url: |-
  https://cgc.sbgenomics.com/ns/brood/images/admin/sbg-public-data/broad-best-practice-data-pre-processing-workflow-4-1-0-0/26.png
sbg:latestRevision: 26
sbg:license: BSD 3-Clause License
sbg:links:
- id: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
  label: Homepage
- id: https://github.com/gatk-workflows/gatk4-data-processing
  label: Source Code
- id: |-
    https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
  label: Download
- id: https://www.ncbi.nlm.nih.gov/pubmed?term=20644199
  label: Publications
- id: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/
  label: Documentation
sbg:modifiedBy: admin
sbg:modifiedOn: 1612280619
sbg:project: admin/sbg-public-data
sbg:projectName: SBG Public data
sbg:publisher: sbg
sbg:revision: 26
sbg:revisionNotes: NO parameters exposed
sbg:revisionsInfo:
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002743
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002743
  sbg:revision: 1
  sbg:revisionNotes: 'dev - v2: labels added, description missing'
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002743
  sbg:revision: 2
  sbg:revisionNotes: v17 - dev project
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002743
  sbg:revision: 3
  sbg:revisionNotes: v18 - dev
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002743
  sbg:revision: 4
  sbg:revisionNotes: Mark Duplicates updated
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002744
  sbg:revision: 5
  sbg:revisionNotes: GatherBamFiles - exposed out_prefix
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002744
  sbg:revision: 6
  sbg:revisionNotes: Add BWA BAM output
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002744
  sbg:revision: 7
  sbg:revisionNotes: Expose smart pairing output in BWA
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002744
  sbg:revision: 8
  sbg:revisionNotes: Revert back to rev5
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002744
  sbg:revision: 9
  sbg:revisionNotes: Expose bwa bam filename
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002744
  sbg:revision: 10
  sbg:revisionNotes: Revert back to rev5
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002745
  sbg:revision: 11
  sbg:revisionNotes: dev - v26
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002745
  sbg:revision: 12
  sbg:revisionNotes: Documentation improved by Marko Marinkovic
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1572002745
  sbg:revision: 13
  sbg:revisionNotes: perf bench updated
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1581357120
  sbg:revision: 14
  sbg:revisionNotes: requrements added - to enable protability
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1581357121
  sbg:revision: 15
  sbg:revisionNotes: dev41
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1581357122
  sbg:revision: 16
  sbg:revisionNotes: dev - v42
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1581357122
  sbg:revision: 17
  sbg:revisionNotes: Fix PL RG issue
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1581360365
  sbg:revision: 18
  sbg:revisionNotes: Remove the default RG PL bit, and add proper description.
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1581524007
  sbg:revision: 19
  sbg:revisionNotes: Remove ignore_rg_information parameter; add requiremnets for
    cwl-tool
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1589907021
  sbg:revision: 20
  sbg:revisionNotes: use soft clipping
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1612280618
  sbg:revision: 21
  sbg:revisionNotes: v60 dev - bwa mem and samtools view
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1612280619
  sbg:revision: 22
  sbg:revisionNotes: v61 - dev
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1612280619
  sbg:revision: 23
  sbg:revisionNotes: v63
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1612280619
  sbg:revision: 24
  sbg:revisionNotes: v65
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1612280619
  sbg:revision: 25
  sbg:revisionNotes: no header
- sbg:modifiedBy: admin
  sbg:modifiedOn: 1612280619
  sbg:revision: 26
  sbg:revisionNotes: NO parameters exposed
sbg:sbgMaintained: false
sbg:toolAuthor: BROAD
sbg:validationErrors: []
sbg:wrapperAuthor: Seven Bridges
