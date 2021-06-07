$namespaces: {dct: 'http://purl.org/dc/terms/', foaf: 'http://xmlns.com/foaf/0.1/',
  sbg: 'https://sevenbridges.com'}
class: Workflow
cwlVersion: v1.0
dct:creator: {'@id': sbg, 'foaf:mbox': 'mailto:support@sbgenomics.com', 'foaf:name': SevenBridges}
doc: "**BROAD Best Practice Data Pre-processing Workflow 4.1.0.0**  is used to prepare
  data for variant calling analysis. \n\nIt can be divided into two major segments:
  alignment to reference genome and data cleanup operations that correct technical
  biases [1].\n\n*A list of all inputs and parameters with corresponding descriptions
  can be found at the bottom of this page.*\n\n### Common Use Cases\n\n* **BROAD Best
  Practice Data Pre-processing Workflow 4.1.0.0**  is designed to operate on individual
  samples.\n* Resulting BAM files are ready for variant calling analysis and can be
  further processed by other BROAD best practice pipelines, like **Generic germline
  short variant per-sample calling workflow** [2], **Somatic CNVs workflow** [3] and
  **Somatic SNVs+Indel workflow** [4].\n\n\n### Changes Introduced by Seven Bridges\n\nThis
  pipeline represents the CWL implementation of BROADs [original WDL file](https://github.com/gatk-workflows/gatk4-data-processing)
  available on github. Minor differences are introduced in order to successfully adapt
  to the Seven Bridges Platform. These differences are listed below:\n* **SamToFastqAndBwaMem**
  step is divided into elementary steps: **SamToFastq** and  **BWA Mem**  \n* **SortAndFixTags**
  is divided into elementary steps: **SortSam** and **SetNmMdAndUqTags**\n* Added
  **SBG Lines to Interval List**: this tool is used to adapt results obtained with
  **CreateSequenceGroupingTSV**  for platform execution, more precisely for scattering.\n\n\n###
  Common Issues and Important Notes\n\n* **BROAD Best Practice Data Pre-processing
  Workflow 4.1.0.0**  expects unmapped BAM file format as the main input.\n* **Input
  Alignments** (`--in_alignments`) - provided unmapped BAM (uBAM) file should be in
  query-sorter order and all reads must have RG tags. Also, input uBAM files must
  pass validation by **ValidateSamFile**.\n* For each tool in the workflow, equivalent
  parameter settings to the one listed in the corresponding WDL file are set as defaults.
  \n\n### Performance Benchmarking\nSince this CWL implementation is meant to be equivalent
  to GATKs original WDL, there are no additional optimization steps beside instance
  and storage definition. \nThe c5.9xlarge AWS instance hint is used for WGS inputs
  and attached storage is set to 1.5TB.\nIn the table given below one can find results
  of test runs for WGS and WES samples. All calculations are performed with reference
  files corresponding to assembly 38.\n\n*Cost can be significantly reduced by spot
  instance usage. Visit knowledge center for more details.*\n\n| Input Size | Experimental
  Strategy | Coverage| Duration | Cost (spot) | AWS Instance Type |\n| --- | --- |
  --- | --- | --- | --- | \n| 6.6 GiB | WES | 70 |1h 9min | $2.02 | c5.9 |\n|3.4 GiB
  | WES |  40 | 39min   | $1.15 | c5.9 |\n|2.1 GiB | WES |  20 | 28min   | $0.82 |
  c5.9 |\n|0.7 GiB | WES |  10 | 14min   | $0.42 | c5.9 |\n| 116 GiB   | WGS (HG001)
  | 50 | 1day 1h 22min   | $49.19 | r4.8 |\n| 185 GiB   | WGS | 50 |1day 8h  14 min
  \  | $56.05 | c5.9 |\n| 111.3 GiB| WGS | 30 |19h 27min | $33.82 | c5.9 |\n| 37.2
  GiB  | WGS | 10 |6h 22min   | $11.09 | c5.9 |\n\n\n\n### API Python Implementation\nThe
  app's draft task can also be submitted via the **API**. In order to learn how to
  get your **Authentication token** and **API endpoint** for corresponding platform
  visit our [documentation](https://github.com/sbg/sevenbridges-python#authentication-and-configuration).\n\n```python\n#
  Initialize the SBG Python API\nfrom sevenbridges import Api\napi = Api(token=\"enter_your_token\",
  url=\"enter_api_endpoint\")\n# Get project_id/app_id from your address bar. Example:
  https://igor.sbgenomics.com/u/your_username/project/app\nproject_id = \"your_username/project\"\napp_id
  = \"your_username/project/app\"\n# Replace inputs with appropriate values\ninputs
  = {\n\t\"in_alignments\": list(api.files.query(project=project_id, names=[\"HCC1143BL.reverted.bam\"])),
  \n\t\"reference_index_tar\": api.files.query(project=project_id, names=[\"Homo_sapiens_assembly38.fasta.tar\"])[0],
  \n\t\"in_reference\": api.files.query(project=project_id, names=[\"Homo_sapiens_assembly38.fasta\"])[0],
  \n\t\"ref_dict\": api.files.query(project=project_id, names=[\"Homo_sapiens_assembly38.dict\"])[0],\n\t\"known_snps\":
  api.files.query(project=project_id, names=[\"1000G_phase1.snps.high_confidence.hg38.vcf\"])[0],\n
  \       \"known_indels\": list(api.files.query(project=project_id, names=[\"Homo_sapiens_assembly38.known_indels.vcf\",
  Mills_and_1000G_gold_standard.indels.hg38.vcf]))}\n# Creates draft task\ntask =
  api.tasks.create(name=\"BROAD Best Practice Data Pre-processing Workflow 4.1.0.0
  - API Run\", project=project_id, app=app_id, inputs=inputs, run=False)\n```\n\nInstructions
  for installing and configuring the API Python client, are provided on [github](https://github.com/sbg/sevenbridges-python#installation).
  For more information about using the API Python client, consult [the client documentation](http://sevenbridges-python.readthedocs.io/en/latest/).
  **More examples** are available [here](https://github.com/sbg/okAPI).\n\nAdditionally,
  [API R](https://github.com/sbg/sevenbridges-r) and [API Java](https://github.com/sbg/sevenbridges-java)
  clients are available. To learn more about using these API clients please refer
  to the [API R client documentation](https://sbg.github.io/sevenbridges-r/), and
  [API Java client documentation](https://docs.sevenbridges.com/docs/java-library-quickstart).\n\n\n###
  References\n\n[1] [Data Pre-processing](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165)\n\n[2]
  [Generic germline short variant per-sample calling](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145)\n\n[3]
  [Somatic CNVs](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11147)\n\n[4]
  [Somatic SNVs+Indel pipeline ](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146)"
hints:
- {class: 'sbg:AWSInstanceType', value: c5.9xlarge;ebs-gp2;3000}
id: uros_sipetic/broad-best-practice-data-pre-processing-4-1-0-0-demo/broad-best-practice-data-pre-processing-workflow-4-1-0-0/14
inputs:
- doc: Input alignments files in unmapped BAM format.
  id: in_alignments
  label: Input alignments
  sbg:fileTypes: SAM, BAM
  sbg:x: -623.1131591796875
  sbg:y: 39.18902587890625
  type: {items: File, type: array}
- {doc: FASTA reference or BWA index archive., id: reference_index_tar, label: BWA
    index archive or FASTA reference, 'sbg:fileTypes': 'TAR, FA, FASTA', 'sbg:x': -612.580078125,
  'sbg:y': 425.3470764160156, type: File}
- doc: Number of threads.
  id: threads
  label: Threads
  sbg:exposed: true
  type: ['null', int]
- doc: Input reference in FASTA format.
  id: in_reference
  label: FASTA reference
  sbg:fileTypes: FASTA, FA
  sbg:x: -612.513916015625
  sbg:y: 555.003173828125
  secondaryFiles: [.fai, ^.dict]
  type: File
- doc: VCF file with known SNPs.
  id: known_snps
  label: Known SNPs
  sbg:fileTypes: VCF
  sbg:x: 618.927734375
  sbg:y: 590.1701049804688
  secondaryFiles: [.idx]
  type:
  - 'null'
  - {items: File, type: array}
- doc: VCF file(s) with known INDELs.
  id: known_indels
  label: Known INDELs
  sbg:fileTypes: VCF
  sbg:x: 617.8395385742188
  sbg:y: 456.47503662109375
  secondaryFiles: [.idx]
  type:
  - 'null'
  - {items: File, type: array}
- {doc: DICT file corresponding to the FASTA reference., id: ref_dict, label: DICT
    file, 'sbg:fileTypes': DICT, 'sbg:x': 599.5844116210938, 'sbg:y': -34.96286392211914,
  type: File}
- doc: When writing files that need to be sorted, this will specify the number of
    records stored in RAM before spilling to disk. Increasing this number reduces
    the number of file handles needed to sort the file, and increases the amount of
    RAM needed.
  id: max_records_in_ram
  label: Max records in RAM
  sbg:exposed: true
  type: ['null', int]
- id: output_prefix
  sbg:exposed: true
  type: ['null', string]
label: BROAD Best Practice Data Pre-processing Workflow 4.1.0.0
outputs:
- doc: Output BAM file.
  id: out_alignments
  label: Output BAM file
  outputSource: [gatk_gatherbamfiles_4_1_0_0/out_alignments]
  sbg:fileTypes: BAM
  sbg:x: 2052.86767578125
  sbg:y: 289.4576416015625
  type: ['null', File]
- doc: MD5 sum of the output BAM file.
  id: out_md5
  label: MD5 file
  outputSource: [gatk_gatherbamfiles_4_1_0_0/out_md5]
  sbg:fileTypes: MD5
  sbg:x: 2048
  sbg:y: 114.24113464355469
  type: ['null', File]
- id: output_metrics
  outputSource: [gatk_markduplicates_4_1_0_0/output_metrics]
  sbg:fileTypes: METRICS
  sbg:x: 457.1893615722656
  sbg:y: -51.47343826293945
  type: File
requirements:
- {class: StepInputExpressionRequirement}
- {class: InlineJavascriptRequirement}
- {class: ScatterFeatureRequirement}
- {class: SubworkflowFeatureRequirement}
sbg:appVersion: [v1.0]
sbg:categories: [Genomics, Alignment]
sbg:content_hash: a323e872efffe3bf7231ad471fa6889eb779314f09bd28bc315b2aa87746b3bf1
sbg:contributors: [uros_sipetic, nens]
sbg:createdBy: nens
sbg:createdOn: 1561370133
sbg:id: uros_sipetic/broad-best-practice-data-pre-processing-4-1-0-0-demo/broad-best-practice-data-pre-processing-workflow-4-1-0-0/14
sbg:image_url: https://igor.sbgenomics.com/ns/brood/images/uros_sipetic/broad-best-practice-data-pre-processing-4-1-0-0-demo/broad-best-practice-data-pre-processing-workflow-4-1-0-0/14.png
sbg:latestRevision: 14
sbg:license: BSD 3-Clause License
sbg:links:
- {id: 'https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165',
  label: Homepage}
- {id: 'https://github.com/gatk-workflows/gatk4-data-processing', label: Source Code}
- {id: 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip',
  label: Download}
- {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
- {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/current/',
  label: Documentation}
sbg:modifiedBy: nens
sbg:modifiedOn: 1572520572
sbg:project: uros_sipetic/broad-best-practice-data-pre-processing-4-1-0-0-demo
sbg:projectName: BROAD Best Practice - Data pre-processing 4.1.0.0 - Demo
sbg:publisher: sbg
sbg:revision: 14
sbg:revisionNotes: requrements added - to enable protability
sbg:revisionsInfo:
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1561370133, 'sbg:revision': 0, 'sbg:revisionNotes': null}
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1561452349, 'sbg:revision': 1, 'sbg:revisionNotes': 'dev
    - v2: labels added, description missing'}
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1563277221, 'sbg:revision': 2, 'sbg:revisionNotes': v17
    - dev project}
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1563282140, 'sbg:revision': 3, 'sbg:revisionNotes': v18
    - dev}
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1566896948, 'sbg:revision': 4, 'sbg:revisionNotes': Mark
    Duplicates updated}
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1566898662, 'sbg:revision': 5, 'sbg:revisionNotes': GatherBamFiles
    - exposed out_prefix}
- {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1567508021, 'sbg:revision': 6,
  'sbg:revisionNotes': Add BWA BAM output}
- {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1567508199, 'sbg:revision': 7,
  'sbg:revisionNotes': Expose smart pairing output in BWA}
- {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1567508283, 'sbg:revision': 8,
  'sbg:revisionNotes': Revert back to rev5}
- {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1567869452, 'sbg:revision': 9,
  'sbg:revisionNotes': Expose bwa bam filename}
- {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1567869845, 'sbg:revision': 10,
  'sbg:revisionNotes': Revert back to rev5}
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1570697472, 'sbg:revision': 11, 'sbg:revisionNotes': dev
    - v26}
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1570800993, 'sbg:revision': 12, 'sbg:revisionNotes': Documentation
    improved by Marko Marinkovic}
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1570802986, 'sbg:revision': 13, 'sbg:revisionNotes': perf
    bench updated}
- {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1572520572, 'sbg:revision': 14, 'sbg:revisionNotes': requrements
    added - to enable protability}
sbg:sbgMaintained: false
sbg:toolAuthor: BROAD
sbg:validationErrors: []
sbg:wrapperAuthor: Seven Bridges
steps:
- id: bwa_mem_bundle_0_7_15
  in:
  - {default: '3', id: verbose_level}
  - {default: true, id: smart_pairing_in_input_fastq}
  - id: input_reads
    source: [gatk_samtofastq_4_1_0_0/out_reads]
  - {default: 100000000, id: num_input_bases_in_each_batch}
  - {default: None, id: deduplication}
  - {default: 16, id: threads, source: threads}
  - {id: reference_index_tar, source: reference_index_tar}
  - {default: BAM, id: output_format}
  - {default: false, id: mapQ_of_suplementary}
  - {default: true, id: ignore_default_rg_id}
  label: BWA MEM Bundle
  out:
  - {id: aligned_reads}
  - {id: dups_metrics}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, prefix: '', shellQuote: false, valueFrom: "${\n    var cmd = \"/bin/bash
        -c \\\"\";\n    return cmd + \" export REF_CACHE=${PWD} ; \";\n}"}
    - {position: 1, prefix: '', shellQuote: false, valueFrom: "${\n    var reference_file
        = inputs.reference_index_tar.path.split('/')[inputs.reference_index_tar.path.split('/').length
        - 1];\n    return 'tar -tvf ' + reference_file + ' 1>&2; tar -xf ' + reference_file
        + ' ; ';\n\n}"}
    - {position: 2, shellQuote: false, valueFrom: bwa}
    - {position: 3, shellQuote: false, valueFrom: mem}
    - {position: 116, prefix: '', separate: false, shellQuote: false, valueFrom: "${\n
        \   ///////////////////////////////////////////\n    ///  BIOBAMBAM BAMSORMADUP
        \  //////////////////////\n    ///////////////////////////////////////////\n\n
        \   function common_substring(a, b) {\n        var i = 0;\n        while (a[i]
        === b[i] && i < a.length) {\n            i = i + 1;\n        }\n\n        return
        a.slice(0, i);\n    }\n\n    // Set output file name\n\n  var input_1 = '';
        \n  var input_2 = ''; \n  var full_name = ''; \n  var namet = ''; \n  var
        reads_size = 0; // Not used because of situations when size does not exist!\n
        \ var GB_1 = 1024 * 1024 * 1024;\n  var suggested_memory = ''; \n  var suggested_cpus
        = ''; \n  var dedup = ''; \n  var sort_path = ''; \n  var indexfilename =
        ''; \n  var out_format = ''; \n  var ref_name_arr = '';  \n  var indexfilename
        = ''; \n  var extension = ''; \n  \n  \n  \n   if (inputs.input_reads[0] instanceof
        Array) {\n        input_1 = inputs.input_reads[0][0]; // scatter mode\n        input_2
        = inputs.input_reads[0][1];\n    } else if (inputs.input_reads instanceof
        Array) {\n        input_1 = inputs.input_reads[0];\n        input_2 = inputs.input_reads[1];\n
        \   } else {\n        input_1 = [].concat(inputs.input_reads)[0];\n        input_2
        = input_1;\n    }\n    \n  full_name = input_1.path.split('/')[input_1.path.split('/').length
        - 1];\n\n    if (inputs.output_name) {\n        namet = inputs.output_name;\n
        \   } else if (inputs.input_reads.length == 1) {\n        namet = full_name;\n
        \       if (namet.slice(-3, namet.length)==='.gz' || namet.slice(-3, namet.length)
        === '.GZ')\n            namet = namet.slice(0, -3);\n        if (namet.slice(-3,
        namet.length)==='.fq' || namet.slice(-3, namet.length) === '.FQ')\n            namet
        = namet.slice(0, -3);\n        if (namet.slice(-6, namet.length)==='.fastq'
        || namet.slice(-6, namet.length) === '.FASTQ')\n            namet = namet.slice(0,
        -6);\n\n    } else {\n        full_name = input_2.path.split('/')[input_2.path.split('/').length
        - 1];\n        namet = common_substring(full_name, full_name);\n\n        if
        (namet.slice(-1, namet.length) === '_' || namet.slice(-1, namet.length) ===
        '.')\n            namet = namet.slice(0, -1);\n        if (namet.slice(-2,
        namet.length) === 'p_' || namet.slice(-1, namet.length) === 'p.')\n            namet
        = namet.slice(0, -2);\n        if (namet.slice(-2, namet.length) === 'P_'
        || namet.slice(-1, namet.length) === 'P.')\n            namet = namet.slice(0,
        -2);\n        if (namet.slice(-3, namet.length) === '_p_' || namet.slice(-3,
        namet.length) === '.p.')\n            namet = namet.slice(0, -3);\n        if
        (namet.slice(-3, namet.length) === '_pe' || namet.slice(-3, namet.length)
        === '.pe')\n            namet = namet.slice(0, -3);\n        if (namet.slice(-2,
        namet.length) === '_R' || namet.slice(-1, namet.length) === '.R')\n            namet
        = namet.slice(0, -2);\n    }\n\n    //////////////////////////\n    // Set
        sort memory size\n\n    if (reads_size < GB_1) {\n        suggested_memory
        = 4;\n        suggested_cpus = 1;\n    } else if (reads_size < 10 * GB_1)
        {\n        suggested_memory = 15;\n        suggested_cpus = 8;\n    } else
        {\n        suggested_memory = 58;\n        suggested_cpus = 31;\n    }\n\n
        \   var total_memory = ''; \n    var sorter_memory = ''; \n    var sorter_memory_string
        = ''; \n    var threads = ''; \n    var MAX_THREADS = ''; \n    var ref_name_arr
        = '';\n    var ref_name = ''; \n        \n    if (!inputs.total_memory) {\n
        \       total_memory = suggested_memory;\n    } else {\n        total_memory
        = inputs.total_memory;\n    }\n\n    // TODO:Rough estimation, should be fine-tuned!\n
        \   if (total_memory > 16) {\n        sorter_memory = parseInt(total_memory
        / 3);\n    } else {\n        sorter_memory = 5;\n    }\n\n    if (inputs.sort_memory)
        {\n        sorter_memory_string = inputs.sort_memory + 'GiB';\n    } else
        \n        sorter_memory_string = sorter_memory + 'GiB';\n\n    // Read number
        of threads if defined\n    if (inputs.sambamba_threads) {\n        threads
        = inputs.sambamba_threads;\n    } else if (inputs.threads) {\n        threads
        = inputs.threads;\n    } else if (inputs.wgs_hg38_mode_threads) {\n        MAX_THREADS
        = 36;\n        ref_name_arr = inputs.reference_index_tar.path.split('/');\n
        \       ref_name = ref_name_arr[ref_name_arr.length - 1];\n        if (ref_name.search('38')
        >= 0) {\n            threads = inputs.wgs_hg38_mode_threads;\n        } else
        {\n            threads = MAX_THREADS;\n        }\n    } else {\n        threads
        = 8;\n    }\n\n    \n        \n        \n  \n\n    if (inputs.deduplication
        == \"MarkDuplicates\") {\n        dedup = ' markduplicates=1';\n    } else
        if (inputs.deduplication == \"RemoveDuplicates\") {\n        dedup = ' rmdup=1';\n
        \   } \n\n    sort_path = 'bamsormadup';\n    \n    // Coordinate Sorted BAM
        is default\n    if (inputs.output_format == 'CRAM') {\n        out_format
        = ' outputformat=cram SO=coordinate';\n        ref_name_arr = inputs.reference_index_tar.path.split('/');\n
        \       ref_name = ref_name_arr[ref_name_arr.length - 1].split('.tar')[0];
        \       \n        out_format += ' reference=' + ref_name;\n        indexfilename
        = ' indexfilename=' + namet + '.cram.crai';\n        extension = '.cram';\n
        \     \n    } else if (inputs.output_format == 'SAM') {\n        out_format
        = ' outputformat=sam SO=coordinate';\n        extension = '.sam';\n      \n
        \   } else if (inputs.output_format == 'Queryname Sorted BAM') {\n        out_format
        = ' outputformat=bam SO=queryname';\n        extension = '.bam';\n      \n
        \   } else if (inputs.output_format == 'Queryname Sorted SAM') {\n        out_format
        = ' outputformat=sam SO=queryname';\n        extension = '.sam';\n      \n
        \   } else {\n        out_format = ' outputformat=bam SO=coordinate';\n        indexfilename
        = ' indexfilename=' + namet + '.bam.bai';\n        extension = '.bam';\n    }\n
        \   var cmd = \" | \" + sort_path + \" threads=8 level=1 tmplevel=-1 inputformat=sam\";\n
        \   cmd += out_format;\n    cmd += indexfilename;\n    // capture metrics
        file\n    cmd += \" M=\" + namet + \".sormadup_metrics.log\";\n    return
        cmd + ' > ' + namet + extension;\n\n}"}
    - {position: 5, prefix: '', shellQuote: false, valueFrom: "${ function add_param(key,
        val) {\n        if (!val) {\n            return '';\n        }\n        param_list.push(key
        + ':' + val);\n    }\n  \n    var param_list = [];\n    var input_1 = '';
        \n    var read_metadata = ''; \n    var rg_platform  = ''; \n    \n    if
        (inputs.read_group_header) {\n        return '-R ' + inputs.read_group_header;\n
        \   }\n    \n\n    // Set output file name\n    if (inputs.input_reads[0]
        instanceof Array) {\n        input_1 = inputs.input_reads[0][0]; // scatter
        mode\n    } else if (inputs.input_reads instanceof Array) {\n        input_1
        = inputs.input_reads[0];\n    } else {\n        input_1 = [].concat(inputs.input_reads)[0];\n
        \   }\n\n    //Read metadata for input reads\n    read_metadata = input_1.metadata;\n
        \ \n    if (!read_metadata) \n      read_metadata = [];\n\n    if (inputs.rg_id)
        {\n        add_param('ID', inputs.rg_id);\n    } else {\n        add_param('ID',
        '1');\n    }\n\n\n    if (inputs.rg_data_submitting_center) {\n        add_param('CN',
        inputs.rg_data_submitting_center);\n    } else if ('data_submitting_center'
        in read_metadata) {\n        add_param('CN', read_metadata.data_submitting_center);\n
        \   }\n\n    if (inputs.rg_library_id) {\n        add_param('LB', inputs.rg_library_id);\n
        \   } else if ('library_id' in read_metadata) {\n        add_param('LB', read_metadata.library_id);\n
        \   }\n\n    if (inputs.rg_median_fragment_length) {\n        add_param('PI',
        inputs.rg_median_fragment_length);\n    }\n\n\n    if (inputs.rg_platform)
        {\n        add_param('PL', inputs.rg_platform);\n    } else if ('platform'
        in read_metadata) {\n        if (read_metadata.platform == 'HiSeq X Ten')
        {\n            rg_platform = 'Illumina';\n        } else {\n            rg_platform
        = read_metadata.platform;\n        }\n        add_param('PL', rg_platform);\n
        \   }\n\n    if (inputs.rg_platform_unit_id) {\n        add_param('PU', inputs.rg_platform_unit_id);\n
        \   } else if ('platform_unit_id' in read_metadata) {\n        add_param('PU',
        read_metadata.platform_unit_id);\n    }\n\n    if (inputs.rg_sample_id) {\n
        \       add_param('SM', inputs.rg_sample_id);\n    } else if ('sample_id'
        in read_metadata) {\n        add_param('SM', read_metadata.sample_id);\n    }\n
        \   \n    if (!inputs.ignore_default_rg_id) {\n      return \"-R '@RG\\\\t\"
        + param_list.join('\\\\t') + \"'\";\n    } else {\n      return '';\n    }\n\n}"}
    - {position: 105, prefix: '', shellQuote: false, valueFrom: "${\n    /////// Set
        input reads in the correct order depending of the paired end from metadata\n\n
        \   // Set output file name\n  \n    var input_reads = ''; \n    var read_metadata
        = ''; \n    var order = 0 ;// Consider this as normal order given at input:
        pe1 pe2\n    var pe1 = ''; \n  \n    if (inputs.input_reads[0] instanceof
        Array) {\n        input_reads = inputs.input_reads[0]; // scatter mode\n    }
        else {\n        input_reads = inputs.input_reads = [].concat(inputs.input_reads);\n
        \   }\n\n    //Read metadata for input reads\n    read_metadata = input_reads[0].metadata;\n
        \   if (!read_metadata) \n      read_metadata = [];    \n    // Check if paired
        end 1 corresponds to the first given read\n    if (read_metadata == []) {\n
        \       order = 0;\n    } else if ('paired_end' in read_metadata) {\n        pe1
        = read_metadata.paired_end;\n        if (pe1 != 1) \n          order = 1;
        // change order\n    }\n\n    // Return reads in the correct order\n    if
        (input_reads.length == 1) {\n        return input_reads[0].path; // Only one
        read present\n    } else if (input_reads.length == 2) {\n        if (order
        == 0) \n          return input_reads[0].path + ' ' + input_reads[1].path;\n
        \       else \n          return input_reads[1].path + ' ' + input_reads[0].path;\n
        \   }\n    return ''; \n\n}"}
    - {position: 6, prefix: -t, shellQuote: false, valueFrom: "${\n    var MAX_THREADS
        = 36;\n    var suggested_threads = 8;\n    var threads = ''; \n    var ref_name_arr
        = '';\n    var ref_name = '';\n\n    if (inputs.threads) {\n        threads
        = inputs.threads;\n    } else if (inputs.wgs_hg38_mode_threads) {\n        ref_name_arr
        = inputs.reference_index_tar.path.split('/');\n        ref_name = ref_name_arr[ref_name_arr.length
        - 1];\n        if (ref_name.search('38') >= 0) {\n            threads = inputs.wgs_hg38_mode_threads;\n
        \       } else {\n            threads = MAX_THREADS;\n        }\n    } else
        {\n        threads = suggested_threads;\n    }\n\n    return threads;\n}"}
    - {position: 14, prefix: '', shellQuote: false, valueFrom: "${\n    var namet
        = '';\n    var reference_file = ''; \n    var metadata = [].concat(inputs.reference_index_tar)[0].metadata;\n\n
        \   if (metadata && metadata.reference_genome) {\n        namet = metadata.reference_genome;\n
        \   } else {\n        reference_file = inputs.reference_index_tar.path.split('/')[inputs.reference_index_tar.path.split('/').length
        - 1];\n        namet = reference_file.slice(0, -4); // cut .tar extension;
        \n    }\n\n    return namet;\n}"}
    - {position: 10004, prefix: '', shellQuote: false, valueFrom: "${\n    var cmd
        = \";declare -i pipe_statuses=(\\\\${PIPESTATUS[*]});len=\\\\${#pipe_statuses[@]};declare
        -i tot=0;echo \\\\${pipe_statuses[*]};for (( i=0; i<\\\\${len}; i++ ));do
        if [ \\\\${pipe_statuses[\\\\$i]} -ne 0 ];then tot=\\\\${pipe_statuses[\\\\$i]};
        fi;done;if [ \\\\$tot -ne 0 ]; then >&2 echo Error in piping. Pipe statuses:
        \\\\${pipe_statuses[*]};fi; if [ \\\\$tot -ne 0 ]; then false;fi\\\"\";\n
        \   return cmd;\n}"}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "**BWA MEM** is an algorithm designed for aligning sequence reads onto a
      large reference genome. BWA MEM is implemented as a component of BWA. The algorithm
      can automatically choose between performing end-to-end and local alignments.
      BWA MEM is capable of outputting multiple alignments, and finding chimeric reads.
      It can be applied to a wide range of read lengths, from 70 bp to several megabases.
      \n\nIn order to obtain possibilities for additional fast processing of aligned
      reads, two tools are embedded together into the same package with BWA MEM (0.7.13):
      Samblaster. (0.1.22) and Sambamba (v0.6.0). \nIf deduplication of alignments
      is needed, it can be done by setting the parameter 'Duplication'. **Samblaster**
      will be used internally to perform this action.\nBesides the standard BWA MEM
      SAM output file, BWA MEM package has been extended to support two additional
      output options: a BAM file obtained by piping through **Sambamba view** while
      filtering out the secondary alignments, as well as a Coordinate Sorted BAM option
      that additionally pipes the output through **Sambamba sort**, along with an
      accompanying .bai file produced by **Sambamba sort** as side effect. Sorted
      BAM is the default output of BWA MEM. Parameters responsible for these additional
      features are 'Filter out secondary alignments' and 'Output format'. Passing
      data from BWA MEM to Samblaster and Sambamba tools has been done through the
      pipes which saves processing times of two read and write of aligned reads into
      the hard drive. \n\nFor input reads fastq files of total size less than 10 GB
      we suggest using the default setting for parameter 'total memory' of 15GB, for
      larger files we suggest using 58 GB of memory and 32 CPU cores.\n\n**Important:**\nIn
      order to work BWA MEM Bundle requires fasta reference file accompanied with
      **bwa fasta indices** in TAR file.\nThere is the **known issue** with samblaster.
      It does not support processing when number of sequences in fasta is larger than
      32768. If this is the case do not use deduplication option because the output
      BAM will be corrupted.\n\nHuman reference genome version 38 comes with ALT contigs,
      a collection of diverged alleles present in some humans but not the others.
      Making effective use of these contigs will help to reduce mapping artifacts,
      however, to facilitate mapping these ALT contigs to the primary assembly, GRC
      decided to add to each contig long flanking sequences almost identical to the
      primary assembly. As a result, a naive mapping against GRCh38+ALT will lead
      to many mapQ-zero mappings in these flanking regions. Please use post-processing
      steps to fix these alignments or implement [steps](https://sourceforge.net/p/bio-bwa/mailman/message/32845712/)
      described by the author of BWA toolkit."
    id: nens/bwa-0-7-15-cwl1-0-demo/bwa-mem-bundle-0-7-15/13
    inputs:
    - doc: Drop chains shorter than FLOAT fraction of the longest overlapping chain.
      id: drop_chains_fraction
      inputBinding: {position: 4, prefix: -D, shellQuote: false}
      label: Drop chains fraction
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '0.50'
      type: ['null', float]
    - doc: 'Verbose level: 1=error, 2=warning, 3=message, 4+=debugging.'
      id: verbose_level
      inputBinding: {position: 4, prefix: -v, shellQuote: false}
      label: Verbose level
      sbg:category: BWA Input/output options
      sbg:toolDefaultValue: '3'
      type:
      - 'null'
      - name: verbose_level
        symbols: ['1', '2', '3', '4']
        type: enum
    - doc: Amount of RAM [Gb] to give to the sorting algorithm (if not provided will
        be set to one third of the total memory).
      id: sort_memory
      label: Memory for BAM sorting
      sbg:category: Execution
      type: ['null', int]
    - doc: Lower the number of threads if HG38 reference genome is used.
      id: wgs_hg38_mode_threads
      label: Optimize threads for HG38
      sbg:category: Execution
      sbg:toolDefaultValue: 'False'
      type: ['null', int]
    - doc: Band width for banded alignment.
      id: band_width
      inputBinding: {position: 4, prefix: -w, shellQuote: false}
      label: Band width
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '100'
      type: ['null', int]
    - doc: Smart pairing in input FASTQ file (ignoring in2.fq).
      id: smart_pairing_in_input_fastq
      inputBinding: {position: 4, prefix: -p, shellQuote: false}
      label: Smart pairing in input FASTQ file
      sbg:category: BWA Input/output options
      type: ['null', boolean]
    - doc: Specify the identifier for the sequencing library preparation, which will
        be placed in RG line.
      id: rg_library_id
      label: Library ID
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Inferred from metadata
      type: ['null', string]
    - doc: Perform at most INT rounds of mate rescues for each read.
      id: mate_rescue_rounds
      inputBinding: {position: 4, prefix: -m, shellQuote: false}
      label: Mate rescue rounds
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '50'
      type: ['null', string]
    - doc: Reserved number of threads on the instance used by scheduler.
      id: reserved_threads
      label: Reserved number of threads on the instance
      sbg:category: Configuration
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    - doc: Input sequence reads.
      id: input_reads
      label: Input reads
      sbg:category: Input files
      sbg:fileTypes: FASTQ, FASTQ.GZ, FQ, FQ.GZ
      type: {items: File, type: array}
    - doc: Penalty for an unpaired read pair.
      id: unpaired_read_penalty
      inputBinding: {position: 4, prefix: -U, shellQuote: false}
      label: Unpaired read penalty
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '17'
      type: ['null', int]
    - doc: Penalty for 5'- and 3'-end clipping.
      id: clipping_penalty
      inputBinding: {itemSeparator: ',', position: 4, prefix: -L, separate: false,
        shellQuote: false}
      label: Clipping penalty
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '[5,5]'
      type:
      - 'null'
      - {items: int, type: array}
    - doc: Look for internal seeds inside a seed longer than {-k} * FLOAT.
      id: select_seeds
      inputBinding: {position: 4, prefix: -r, shellQuote: false}
      label: Select seeds
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '1.5'
      type: ['null', float]
    - doc: Score for a sequence match, which scales options -TdBOELU unless overridden.
      id: score_for_a_sequence_match
      inputBinding: {position: 4, prefix: -A, shellQuote: false}
      label: Score for a sequence match
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    - doc: Off-diagonal X-dropoff.
      id: dropoff
      inputBinding: {position: 4, prefix: -d, shellQuote: false}
      label: Dropoff
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '100'
      type: ['null', int]
    - doc: process INT input bases in each batch regardless of nThreads (for reproducibility)
      id: num_input_bases_in_each_batch
      inputBinding: {position: 4, prefix: -K, shellQuote: false}
      label: process INT input bases in each batch (for reproducibility)
      type: ['null', int]
    - doc: Total memory to be used by the tool in GB. It's sum of BWA, Sambamba Sort
        and Samblaster. For fastq files of total size less than 10GB, we suggest using
        the default setting of 15GB, for larger files we suggest using 58GB of memory
        (and 32CPU cores).
      id: total_memory
      label: Total memory
      sbg:category: Execution
      sbg:toolDefaultValue: '15'
      type: ['null', int]
    - doc: Gap extension penalty; a gap of size k cost '{-O} + {-E}*k'. This array
        can't have more than two values.
      id: gap_extension_penalties
      inputBinding: {itemSeparator: ',', position: 4, prefix: -E, separate: false,
        shellQuote: false}
      label: Gap extension
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '[1,1]'
      type:
      - 'null'
      - {items: int, type: array}
    - doc: Use Samblaster for finding duplicates on sequence reads.
      id: deduplication
      label: PCR duplicate detection
      sbg:category: Samblaster parameters
      sbg:toolDefaultValue: MarkDuplicates
      type:
      - 'null'
      - name: deduplication
        symbols: [None, MarkDuplicates, RemoveDuplicates]
        type: enum
    - doc: Treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt
        file).
      id: ignore_alt_file
      inputBinding: {position: 4, prefix: -j, shellQuote: false}
      label: Ignore ALT file
      sbg:category: BWA Input/output options
      type: ['null', boolean]
    - doc: Read group ID
      id: rg_id
      label: Read group ID
      sbg:category: Configuration
      sbg:toolDefaultValue: '1'
      type: ['null', string]
    - doc: Use soft clipping for supplementary alignments.
      id: use_soft_clipping
      inputBinding: {position: 4, prefix: -Y, shellQuote: false}
      label: Use soft clipping
      sbg:category: BWA Input/output options
      type: ['null', boolean]
    - doc: If there are <INT hits with score >80% of the max score, output all in
        XA. This array should have no more than two values.
      id: output_in_xa
      inputBinding: {itemSeparator: ',', position: 4, prefix: -h, separate: false,
        shellQuote: false}
      label: Output in XA
      sbg:category: BWA Input/output options
      sbg:toolDefaultValue: '[5, 200]'
      type:
      - 'null'
      - {items: int, type: array}
    - doc: Specify the version of the technology that was used for sequencing, which
        will be placed in RG line.
      id: rg_platform
      label: Platform
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Inferred from metadata
      type:
      - 'null'
      - name: rg_platform
        symbols: ['454', Helicos, Illumina, Solid, IonTorrent]
        type: enum
    - doc: Number of threads for BWA, Samblaster and Sambamba sort process.
      id: threads
      label: Threads
      sbg:category: Execution
      sbg:toolDefaultValue: '8'
      type: ['null', int]
    - doc: Skip pairing; mate rescue performed unless -S also in use.
      id: skip_pairing
      inputBinding: {position: 4, prefix: -P, shellQuote: false}
      label: Skip pairing
      sbg:category: BWA Algorithm options
      type: ['null', boolean]
    - doc: Insert STR to header if it starts with @; or insert lines in FILE.
      id: insert_string_to_header
      inputBinding: {position: 4, prefix: -H, shellQuote: false}
      label: Insert string to output SAM or BAM header
      sbg:category: BWA Input/output options
      type: ['null', string]
    - doc: Output the reference FASTA header in the XR tag.
      id: output_header
      inputBinding: {position: 4, prefix: -V, shellQuote: false}
      label: Output header
      sbg:category: BWA Input/output options
      type: ['null', boolean]
    - doc: Seed occurrence for the 3rd round seeding.
      id: seed_occurrence_for_the_3rd_round
      inputBinding: {position: 4, prefix: -y, shellQuote: false}
      label: Seed occurrence for the 3rd round
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '20'
      type: ['null', int]
    - doc: 'Sequencing technology-specific settings; Setting -x changes multiple parameters
        unless overriden. pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads
        to ref). ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads
        to ref). intractg: -B9 -O16 -L5  (intra-species contigs to ref).'
      id: read_type
      inputBinding: {position: 4, prefix: -x, shellQuote: false}
      label: Sequencing technology-specific settings
      sbg:category: BWA Scoring options
      type:
      - 'null'
      - name: read_type
        symbols: [pacbio, ont2d, intractg]
        type: enum
    - {doc: Reference fasta file with BWA index files packed in TAR., id: reference_index_tar,
      label: Reference Index TAR, 'sbg:category': Input files, 'sbg:fileTypes': TAR,
      type: File}
    - doc: Mark shorter split hits as secondary.
      id: mark_shorter
      inputBinding: {position: 4, prefix: -M, shellQuote: false}
      label: Mark shorter
      sbg:category: BWA Input/output options
      type: ['null', boolean]
    - doc: Specify the mean, standard deviation (10% of the mean if absent), max (4
        sigma from the mean if absent) and min of the insert size distribution.FR
        orientation only. This array can have maximum four values, where first two
        should be specified as FLOAT and last two as INT.
      id: speficy_distribution_parameters
      inputBinding: {itemSeparator: ' -I', position: 4, prefix: -I, separate: false,
        shellQuote: false}
      label: Specify distribution parameters
      sbg:category: BWA Input/output options
      type:
      - 'null'
      - {items: float, type: array}
    - doc: Minimum alignment score for a read to be output in SAM/BAM.
      id: minimum_output_score
      inputBinding: {position: 4, prefix: -T, shellQuote: false}
      label: Minimum alignment score for a read to be output in SAM/BAM
      sbg:category: BWA Input/output options
      sbg:toolDefaultValue: '30'
      type: ['null', int]
    - doc: Cordinate sort is default output.
      id: output_format
      label: Output format
      sbg:category: Execution
      sbg:toolDefaultValue: Coordinate Sorted BAM
      type:
      - 'null'
      - name: output_format
        symbols: [SAM, BAM, CRAM, Queryname Sorted BAM, Queryname Sorted SAM]
        type: enum
    - doc: Skip mate rescue.
      id: skip_mate_rescue
      inputBinding: {position: 4, prefix: -S, shellQuote: false}
      label: Skip mate rescue
      sbg:category: BWA Algorithm options
      type: ['null', boolean]
    - doc: Skip seeds with more than INT occurrences.
      id: skip_seeds
      inputBinding: {position: 4, prefix: -c, shellQuote: false}
      label: Skip seeds with more than INT occurrences
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '500'
      type: ['null', int]
    - doc: Filter out secondary alignments. Sambamba view tool will be used to perform
        this internally.
      id: filter_out_secondary_alignments
      label: Filter out secondary alignments
      sbg:category: Execution
      sbg:toolDefaultValue: 'False'
      type: ['null', boolean]
    - doc: Name of the output BAM file.
      id: output_name
      label: Output SAM/BAM file name
      sbg:category: Configuration
      type: ['null', string]
    - doc: Minimum seed length for BWA MEM.
      id: minimum_seed_length
      inputBinding: {position: 4, prefix: -k, shellQuote: false}
      label: Minimum seed length
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '19'
      type: ['null', int]
    - doc: Gap open penalties for deletions and insertions. This array can't have
        more than two values.
      id: gap_open_penalties
      inputBinding: {itemSeparator: ',', position: 4, prefix: -O, separate: false,
        shellQuote: false}
      label: Gap open penalties
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '[6,6]'
      type:
      - 'null'
      - {items: int, type: array}
    - doc: Specify the median fragment length for RG line.
      id: rg_median_fragment_length
      label: Median fragment length
      sbg:category: BWA Read Group Options
      type: ['null', string]
    - doc: Penalty for a mismatch.
      id: mismatch_penalty
      inputBinding: {position: 4, prefix: -B, shellQuote: false}
      label: Mismatch penalty
      sbg:category: BWA Scoring options
      sbg:toolDefaultValue: '4'
      type: ['null', int]
    - doc: Output all alignments for SE or unpaired PE.
      id: output_alignments
      inputBinding: {position: 4, prefix: -a, shellQuote: false}
      label: Output alignments
      sbg:category: BWA Input/output options
      type: ['null', boolean]
    - doc: Discard full-length exact matches.
      id: discard_exact_matches
      inputBinding: {position: 4, prefix: -e, shellQuote: false}
      label: Discard exact matches
      sbg:category: BWA Algorithm options
      type: ['null', boolean]
    - doc: Specify the platform unit (lane/slide) for RG line - An identifier for
        lanes (Illumina), or for slides (SOLiD) in the case that a library was split
        and ran over multiple lanes on the flow cell or slides.
      id: rg_platform_unit_id
      label: Platform unit ID
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Inferred from metadata
      type: ['null', string]
    - doc: Don't modify mapQ of supplementary alignments
      id: mapQ_of_suplementary
      inputBinding: {position: 4, prefix: -q, shellQuote: false}
      label: Don't modify mapQ of supplementary alignments
      type: ['null', boolean]
    - doc: Specify the sample ID for RG line - A human readable identifier for a sample
        or specimen, which could contain some metadata information. A sample or specimen
        is material taken from a biological entity for testing, diagnosis, propagation,
        treatment, or research purposes, including but not limited to tissues, body
        fluids, cells, organs, embryos, body excretory products, etc.
      id: rg_sample_id
      label: Sample ID
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Inferred from metadata
      type: ['null', string]
    - doc: Specify the data submitting center for RG line.
      id: rg_data_submitting_center
      label: Data submitting center
      sbg:category: BWA Read Group Options
      type: ['null', string]
    - doc: Discard a chain if seeded bases shorter than INT.
      id: discard_chain_length
      inputBinding: {position: 4, prefix: -W, shellQuote: false}
      label: Discard chain length
      sbg:category: BWA Algorithm options
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: for split alignment, take the alignment with the smallest coordinate as
        primary.
      id: split_alignment_primary
      inputBinding: {position: 4, prefix: '-5', shellQuote: false}
      label: Split alignment smallest coordinate as primary
      type: ['null', boolean]
    - doc: Append FASTA/FASTQ comment to SAM output.
      id: append_comment
      inputBinding: {position: 4, prefix: -C, shellQuote: false}
      label: Append comment
      sbg:category: BWA Input/output options
      type: ['null', boolean]
    - doc: Read group header line such as '@RG\tID:foo\tSM:bar'.  This value takes
        precedence over per-attribute parameters.
      id: read_group_header
      label: Read group header
      sbg:category: BWA Read Group Options
      sbg:toolDefaultValue: Constructed from per-attribute parameters or inferred
        from metadata.
      type: ['null', string]
    - doc: Ignore default RG ID ('1').
      id: ignore_default_rg_id
      label: Ignore default RG ID
      sbg:category: BWA Read Group Options
      type: ['null', boolean]
    label: BWA MEM Bundle
    outputs:
    - doc: Aligned reads.
      id: aligned_reads
      label: Aligned SAM/BAM
      outputBinding: {glob: '*am', outputEval: "${  \n    self = inheritMetadata(self,
          inputs.input_reads);\n    for (var i = 0; i < self.length; i++) {\n        var
          met = inputs.reference_index_tar.path.split('/').pop().replace(/.tar/i,'').replace(/.gz/i,'');
          \ // cut .tar and gz extension \n        \n        var out_metadata = {'reference_genome':
          met\n        };\n        \n        self[i] = setMetadata(self[i], out_metadata);\n
          \   };\n\n    return self\n\n}"}
      sbg:fileTypes: SAM, BAM, CRAM
      secondaryFiles: [.bai, ^.bai, .crai, ^.crai]
      type: ['null', File]
    - doc: Metrics file for biobambam mark duplicates
      id: dups_metrics
      label: Sormadup metrics
      outputBinding: {glob: '*.sormadup_metrics.log'}
      sbg:fileTypes: LOG
      type: ['null', File]
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: "${\n    // Calculate suggested number
        of CPUs depending of the input reads size\n    if (inputs.input_reads.constructor
        == Array) {\n        if (inputs.input_reads[1]) {\n            var reads_size
        = inputs.input_reads[0].size + inputs.input_reads[1].size;\n        } else
        {\n            reads_size = inputs.input_reads[0].size;\n        }\n    }
        else {\n        reads_size = inputs.input_reads.size;\n    }\n    if (!reads_size)
        {\n        reads_size = 0;\n    }\n\n\n    var GB_1 = 1024 * 1024 * 1024;\n
        \   if (reads_size < GB_1) {\n        var suggested_cpus = 1;\n    } else
        if (reads_size < 10 * GB_1) {\n        suggested_cpus = 8;\n    } else {\n
        \       suggested_cpus = 31;\n    }\n\n    if (inputs.reserved_threads) {\n
        \       return inputs.reserved_threads;\n    } else if (inputs.threads) {\n
        \       return inputs.threads;\n    } else if (inputs.sambamba_threads) {\n
        \       return inputs.sambamba_threads;\n    } else {\n        return suggested_cpus;\n
        \   }\n}", ramMin: "${\n\n    // Calculate suggested number of CPUs depending
        of the input reads size\n    if (inputs.input_reads.constructor == Array)
        {\n        if (inputs.input_reads[1]) {\n            var reads_size = inputs.input_reads[0].size
        + inputs.input_reads[1].size;\n        } else {\n            reads_size =
        inputs.input_reads[0].size;\n        }\n    } else {\n        reads_size =
        inputs.input_reads.size;\n    }\n    if (!reads_size) {\n        reads_size
        = 0;\n    }\n\n    var GB_1 = 1024 * 1024 * 1024;\n    if (reads_size < GB_1)
        {\n        var suggested_memory = 4;\n    } else if (reads_size < 10 * GB_1)
        {\n        var suggested_memory = 15;\n    } else {\n        var suggested_memory
        = 58;\n    }\n\n    if (inputs.total_memory) {\n        return inputs.total_memory
        * 1024;\n    } else if (inputs.sort_memory) {\n        return inputs.sort_memory
        * 1024;\n    } else {\n        return suggested_memory * 1024;\n    }\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/nens/bwa-0-7-15:0'}
    - class: InitialWorkDirRequirement
      listing: [$(inputs.reference_index_tar), $(inputs.input_reads)]
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n
          \   for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n
          \   }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n
          \   var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2
          = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example
          = o2[i]['metadata'];\n        for (var key in example) {\n            if
          (i == 0)\n                commonMetadata[key] = example[key];\n            else
          {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete
          commonMetadata[key]\n                }\n            }\n        }\n    }\n
          \   if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n
          \   } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i]
          = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar
          toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy
          = function(files, key) {\n    var groupedFiles = [];\n    var tempDict =
          {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n
          \       if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Genomics, Alignment]
    sbg:cmdPreview: '/bin/bash -c " export REF_CACHE=${PWD} ;  tar -tvf reference.HG38.fasta.gz.tar
      1>&2; tar -xf reference.HG38.fasta.gz.tar ;  bwa mem  -R ''@RG\tID:1\tPL:Illumina\tSM:dnk_sample''
      -t 10  reference.HG38.fasta.gz  /path/to/LP6005524-DNA_C01_lane_7.sorted.converted.filtered.pe_2.gz
      /path/to/LP6005524-DNA_C01_lane_7.sorted.converted.filtered.pe_1.gz  | bamsormadup
      threads=8 level=1 tmplevel=-1 inputformat=sam outputformat=cram SO=coordinate
      reference=reference.HG38.fasta.gz indexfilename=LP6005524-DNA_C01_lane_7.sorted.converted.filtered.cram.crai
      M=LP6005524-DNA_C01_lane_7.sorted.converted.filtered.sormadup_metrics.log >
      LP6005524-DNA_C01_lane_7.sorted.converted.filtered.cram  ;declare -i pipe_statuses=(\${PIPESTATUS[*]});len=\${#pipe_statuses[@]};declare
      -i tot=0;echo \${pipe_statuses[*]};for (( i=0; i<\${len}; i++ ));do if [ \${pipe_statuses[\$i]}
      -ne 0 ];then tot=\${pipe_statuses[\$i]}; fi;done;if [ \$tot -ne 0 ]; then >&2
      echo Error in piping. Pipe statuses: \${pipe_statuses[*]};fi; if [ \$tot -ne
      0 ]; then false;fi"'
    sbg:content_hash: adaf2fe4c60a6786bbe84d0896fabd2920d1cac9618c3025407140a9640f5f532
    sbg:contributors: [nens, uros_sipetic]
    sbg:copyOf: nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/22
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1555689212
    sbg:expand_workflow: false
    sbg:id: nens/bwa-0-7-15-cwl1-0-demo/bwa-mem-bundle-0-7-15/13
    sbg:image_url: null
    sbg:latestRevision: 13
    sbg:license: 'BWA: GNU Affero General Public License v3.0, MIT License. Sambamba:
      GNU GENERAL PUBLIC LICENSE. Samblaster: The MIT License (MIT)'
    sbg:links:
    - {id: 'http://bio-bwa.sourceforge.net/', label: Homepage}
    - {id: 'https://github.com/lh3/bwa', label: Source code}
    - {id: 'http://bio-bwa.sourceforge.net/bwa.shtml', label: Wiki}
    - {id: 'http://sourceforge.net/projects/bio-bwa/', label: Download}
    - {id: 'http://arxiv.org/abs/1303.3997', label: Publication}
    - {id: 'http://www.ncbi.nlm.nih.gov/pubmed/19451168', label: Publication BWA Algorithm}
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1558441939
    sbg:project: nens/bwa-0-7-15-cwl1-0-demo
    sbg:projectName: BWA 0.7.15 - cwl1.0- DEMO
    sbg:publisher: sbg
    sbg:revision: 13
    sbg:revisionNotes: Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/22
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1555689212, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/1}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1556035789, 'sbg:revision': 1,
      'sbg:revisionNotes': Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/3}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1556037315, 'sbg:revision': 2,
      'sbg:revisionNotes': Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/4}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1556192655, 'sbg:revision': 3,
      'sbg:revisionNotes': Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/5}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1556193727, 'sbg:revision': 4,
      'sbg:revisionNotes': Copy of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/6}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000453, 'sbg:revision': 5, 'sbg:revisionNotes': Copy
        of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/9}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558002186, 'sbg:revision': 6, 'sbg:revisionNotes': Copy
        of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/10}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558021975, 'sbg:revision': 7, 'sbg:revisionNotes': Copy
        of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/12}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558023132, 'sbg:revision': 8, 'sbg:revisionNotes': Copy
        of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/13}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558085159, 'sbg:revision': 9, 'sbg:revisionNotes': Copy
        of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/15}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558349205, 'sbg:revision': 10, 'sbg:revisionNotes': Copy
        of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/16}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351490, 'sbg:revision': 11, 'sbg:revisionNotes': Copy
        of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/17}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558427784, 'sbg:revision': 12, 'sbg:revisionNotes': Copy
        of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/18}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558441939, 'sbg:revision': 13, 'sbg:revisionNotes': Copy
        of nens/bwa-0-7-15-cwl1-dev/bwa-mem-bundle-0-7-15/22}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Heng Li
    sbg:toolkit: BWA
    sbg:toolkitVersion: 0.7.15
    sbg:validationErrors: []
  sbg:x: -194.6957550048828
  sbg:y: 289.4698181152344
  scatter: [input_reads]
- id: gatk_applybqsr_4_1_0_0
  in:
  - {default: 'true', id: add_output_sam_program_record}
  - {id: bqsr_recal_file, source: gatk_gatherbqsrreports_4_1_0_0/out_gathered_bqsr_reports}
  - {id: in_alignments, source: gatk_setnmmdanduqtags_4_1_0_0/out_alignments}
  - {id: include_intervals_file, source: sbg_lines_to_interval_list_1/out_intervals}
  - {id: in_reference, source: in_reference}
  - default: [10, 20, 30]
    id: static_quantized_quals
  - {default: true, id: use_original_qualities}
  label: GATK ApplyBQSR
  out:
  - {id: out_alignments}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, prefix: '', shellQuote: false, valueFrom: /opt/gatk}
    - {position: 1, shellQuote: false, valueFrom: --java-options}
    - {position: 2, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)
        {\n        return '\\\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\\\"';\n
        \   }\n    return '\\\"-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps
        -XX:+PrintGCDetails -Xloggc:gc_log.log -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10
        -Xms3000m\\\"';\n}"}
    - {position: 3, shellQuote: false, valueFrom: ApplyBQSR}
    - {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n    var in_alignments
        = inputs.in_alignments;\n    var output_ext = inputs.output_file_format ?
        inputs.output_file_format : in_alignments.path.split('.').pop();\n    var
        output_prefix = '';\n    if (inputs.output_prefix)\n    {\n        output_prefix
        = inputs.output_prefix;\n    }\n    else \n    {\n        if (in_alignments.metadata
        && in_alignments.metadata.sample_id)\n        {\n            output_prefix
        = in_alignments.metadata.sample_id;\n        }\n        else \n        {\n
        \           output_prefix = in_alignments.path.split('/').pop().split('.')[0];\n
        \       }\n    }\n    \n    return \"--output \" + output_prefix + \".recalibrated.\"
        + output_ext;\n}"}
    - {position: 5, prefix: '', shellQuote: false, valueFrom: "${\n    var output_ext
        = inputs.output_file_format ? inputs.output_file_format :\n    inputs.in_alignments.path.split('.').pop();\n
        \   if (output_ext == \"cram\")\n    {\n        return '&& for i in *cram.bai;
        do mv \"$i\" \"${i%.cram.bai}.cram.crai\";  done';\n    }\n    else {\n        return
        ''; \n    } \n}"}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "The **GATK ApplyBQSR** tool recalibrates the base qualities score of the
      input reads and outputs a recalibrated BAM, SAM or CRAM file. \n\nThis tool
      performs the second pass in a two-stage process called Base Quality Score Recalibration
      (BQSR). Specifically, it recalibrates the base qualities of the input reads
      based on the recalibration table produced by the **BaseRecalibrator** tool.
      The goal of this procedure is to correct systematic biasses that affect the
      assignment of base quality scores by the sequencer. The first pass consists
      of calculating error empirically and finding patterns in how error varies with
      basecall features over all bases. The relevant observations are written to the
      recalibration table. The second pass consists of applying numerical corrections
      to each individual basecall based on the patterns identified in the first step
      (recorded in the recalibration table) and write out the recalibrated data to
      a new BAM, SAM or CRAM file [1].\n\n*A list of **all inputs and parameters**
      with corresponding descriptions can be found at the bottom of the page.*\n\n###Common
      Use Cases\n\n* The **GATK ApplyBQSR** tool requires the BAM, SAM or CRAM file
      on its **Input BAM/SAM/CRAM file** (`--input`) input and the covariates table
      (= recalibration file) generated by the **BaseRecalibrator** tool on its **BQSR
      recal file** input (`--bqsr-recal-file`). The tool generates on its **Output
      BAM/SAM/CRAM** output a new alignments file which contains recalibrated read
      data.\n\n* Usage example\n\n```\n gatk ApplyBQSR \\\n   --reference reference.fasta
      \\\n   --input input.bam \\\n   --bqsr-recal-file recalibration.table \\\n   --output
      output.bam\n\n```\n\n* If the input alignments file is in CRAM format, the reference
      sequence is required on the **Reference sequence** (`--reference`) input of
      the tool.\n\n* Original qualities can be retained in the output file under the
      \"OQ\" tag if desired. See the **Emit original quals** (`--emit-original-quals`)
      argument for details.\n\n###Changes Introduced by Seven Bridges\n\n* All output
      files will be prefixed using the **Output prefix** parameter. In case **Output
      prefix** is not provided, output prefix will be the same as the Sample ID metadata
      from **Input SAM/BAM/CRAM file**, if the Sample ID metadata exists. Otherwise,
      output prefix will be inferred from the **Input SAM/BAM/CRAM** filename. This
      way, having identical names of the output files between runs is avoided. Moreover,
      \ **recalibrated** will be added before the extension of the output file name.
      \n\n* The user has a possibility to specify the output file format using the
      **Output file format** argument. Otherwise, the output file format will be the
      same as the format of the input file.\n\n* **Include intervals** (`--intervals`)
      option is divided into **Include intervals string** and **Include intervals
      file** options.\n\n* **Exclude intervals** (`--exclude-intervals`) option is
      divided into **Exclude intervals string** and **Exclude intervals file** options.\n\n###Common
      Issues and Important Notes\n\n* Note: This tool replaces the use of PrintReads
      for the application of base quality score recalibration as practiced in earlier
      versions of GATK (2.x and 3.x) [1].\n* Note: You should only run **ApplyBQSR**
      with the covariates table created from the input BAM, SAM or CRAM file [1].\nInput
      BAM/SAM/CRAM file\n* Note: If **Read filter** (`--read-filter`) option is set
      to \"LibraryReadFilter\", **Library** (`--library`) option must be set to some
      value.\n* Note: If **Read filter** (`--read-filter`) option is set to \"PlatformReadFilter\",
      **Platform filter name** (`--platform-filter-name`) option must be set to some
      value.\n* Note: If **Read filter** (`--read-filter`) option is set to\"PlatformUnitReadFilter\",
      **Black listed lanes** (`--black-listed-lanes`) option must be set to some value.
      \n* Note: If **Read filter** (`--read-filter`) option is set to \"ReadGroupBlackListReadFilter\",
      **Read group black list** (`--read-group-black-list`) option must be set to
      some value.\n* Note: If **Read filter** (`--read-filter`) option is set to \"ReadGroupReadFilter\",
      **Keep read group** (`--keep-read-group`) option must be set to some value.\n*
      Note: If **Read filter** (`--read-filter`) option is set to \"ReadLengthReadFilter\",
      **Max read length** (`--max-read-length`) option must be set to some value.\n*
      Note: If **Read filter** (`--read-filter`) option is set to \"ReadNameReadFilter\",
      **Read name** (`--read-name`) option must be set to some value.\n* Note: If
      **Read filter** (`--read-filter`) option is set to \"ReadStrandFilter\", **Keep
      reverse strand only** (`--keep-reverse-strand-only`) option must be set to some
      value.\n* Note: If **Read filter** (`--read-filter`) option is set to \"SampleReadFilter\",
      **Sample** (`--sample`) option must be set to some value.\n\n###Performance
      Benchmarking\n\nBelow is a table describing runtimes and task costs of **GATK
      ApplyBQSR** for a couple of different samples, executed on the AWS cloud instances:\n\n|
      Experiment type |  Input size | Duration |  Cost | Instance (AWS) | \n|:--------------:|:------------:|:--------:|:-------:|:---------:|\n|
      \    RNA-Seq     |  2.2 GB |   7min   | ~0.05$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq
      \    |  6.6 GB |   22min   | ~0.12$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq
      \    | 11 GB |  36min  | ~0.17$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq     |
      22 GB |  1h 15min  | ~0.29$ | c4.2xlarge (8 CPUs) |\n\n*Cost can be significantly
      reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances)
      for more details.*\n\n###References\n\n[1] [GATK ApplyBQSR](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php)"
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-applybqsr-4-1-0-0/10
    inputs:
    - doc: If true, adds a PG tag to created SAM/BAM/CRAM files.
      id: add_output_sam_program_record
      inputBinding: {position: 4, prefix: --add-output-sam-program-record, shellQuote: false}
      label: Add output SAM program record
      sbg:altPrefix: -add-output-sam-program-record
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: add_output_sam_program_record
        symbols: ['true', 'false']
        type: enum
    - doc: Threshold number of ambiguous bases. If null, uses threshold fraction;
        otherwise, overrides threshold fraction. Cannot be used in conjuction with
        argument maxAmbiguousBaseFraction. Valid only if "AmbiguousBaseReadFilter"
        is specified.
      id: ambig_filter_bases
      inputBinding: {position: 4, prefix: --ambig-filter-bases, shellQuote: false}
      label: Ambig filter bases
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', int]
    - doc: Threshold fraction of ambiguous bases. Cannot be used in conjuction with
        argument maxAmbiguousBases. Valid only if "AmbiguousBaseReadFilter" is specified.
      id: ambig_filter_frac
      inputBinding: {position: 4, prefix: --ambig-filter-frac, shellQuote: false}
      label: Ambig filter frac
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', float]
    - doc: Platform unit (PU) to filter out. Valid only if "PlatformUnitReadFilter"
        is specified.
      id: black_listed_lanes
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--black-listed-lanes',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Black listed lanes
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Input recalibration table for BQSR.
      id: bqsr_recal_file
      inputBinding: {position: 4, prefix: --bqsr-recal-file, shellQuote: false}
      label: BQSR recal file
      sbg:altPrefix: -bqsr
      sbg:category: Required Arguments
      sbg:fileTypes: CSV
      type: File
    - doc: If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM
        file.
      id: create_output_bam_index
      inputBinding: {position: 4, prefix: --create-output-bam-index, shellQuote: false}
      label: Create output BAM/CRAM index
      sbg:altPrefix: -OBI
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: create_output_bam_index
        symbols: ['true', 'false']
        type: enum
    - doc: Read filters to be disabled before analysis.
      id: disable_read_filter
      inputBinding: {position: 4, separate: false, shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--disable-read-filter',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Disable read filter
      sbg:altPrefix: -DF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - items:
          name: disable_read_filter
          symbols: [AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter,
            AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter,
            FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter,
            LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter,
            MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter,
            MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter,
            MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter,
            NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter,
            NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter,
            OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter,
            PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter,
            ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter,
            ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter,
            SecondOfPairReadFilter, SeqIsStoredReadFilter, ValidAlignmentEndReadFilter,
            ValidAlignmentStartReadFilter, WellformedReadFilter]
          type: enum
        type: array
    - doc: If specified, do not check the sequence dictionaries from our inputs for
        compatibility. Use at your own risk!
      id: disable_sequence_dictionary_validation
      inputBinding: {position: 4, prefix: --disable-sequence-dictionary-validation,
        shellQuote: false}
      label: Disable sequence dictionary validation
      sbg:altPrefix: -disable-sequence-dictionary-validation
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: 'Disable all tool default read filters (warning: many tools will not function
        correctly without their default read filters on).'
      id: disable_tool_default_read_filters
      inputBinding: {position: 4, prefix: --disable-tool-default-read-filters, shellQuote: false}
      label: Disable tool default read filters
      sbg:altPrefix: -disable-tool-default-read-filters
      sbg:category: Advanced Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Allow a read to be filtered out based on having only 1 soft-clipped block.
        By default, both ends must have a soft-clipped block, setting this flag requires
        only 1 soft-clipped block. Valid only if "OverclippedReadFilter" is specified.
      id: dont_require_soft_clips_both_ends
      inputBinding: {position: 4, prefix: --dont-require-soft-clips-both-ends, shellQuote: false}
      label: Dont require soft clips both ends
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Emit original base qualities under the OQ tag.
      id: emit_original_quals
      inputBinding: {position: 4, prefix: --emit-original-quals, shellQuote: false}
      label: Emit original quals
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: File which contains one or more genomic intervals to exclude from processing.
      id: exclude_intervals_file
      inputBinding: {position: 4, prefix: --exclude-intervals, shellQuote: false}
      label: Exclude intervals file
      sbg:altPrefix: -XL
      sbg:category: Optional Arguments
      sbg:fileTypes: BED
      sbg:toolDefaultValue: 'null'
      type: ['null', File]
    - doc: One or more genomic intervals to exclude from processing.
      id: exclude_intervals_string
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--exclude-intervals',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Exclude intervals string
      sbg:altPrefix: -XL
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Minimum number of aligned bases. Valid only if "OverclippedReadFilter"
        is specified.
      id: filter_too_short
      inputBinding: {position: 4, prefix: --filter-too-short, shellQuote: false}
      label: Filter too short
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '30'
      type: ['null', int]
    - doc: Global Qscore Bayesian prior to use for BQSR.
      id: global_qscore_prior
      inputBinding: {position: 4, prefix: --global-qscore-prior, shellQuote: false}
      label: Global Qscore prior
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '-1.0'
      type: ['null', float]
    - doc: BAM/SAM/CRAM file(s) containing reads.
      id: in_alignments
      inputBinding: {position: 4, prefix: --input, shellQuote: false}
      label: Input BAM/SAM/CRAM file
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM, CRAM
      secondaryFiles: ["${\n    \n    if (self.nameext == \".bam\")\n    {\n        return
          \ self.nameroot + \".bai\";\n    }\n    else if (self.nameext == \".cram\")\n
          \   {\n        return self.nameroot + \".crai\";\n    }\n    \n    return
          ''; \n    \n}"]
      type: File
    - doc: Amount of padding (in bp) to add to each interval you are excluding.
      id: interval_exclusion_padding
      inputBinding: {position: 4, prefix: --interval-exclusion-padding, shellQuote: false}
      label: Interval exclusion padding
      sbg:altPrefix: -ixp
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: Interval merging rule for abutting intervals.
      id: interval_merging_rule
      inputBinding: {position: 4, prefix: --interval-merging-rule, shellQuote: false}
      label: Interval merging rule
      sbg:altPrefix: -imr
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: ALL
      type:
      - 'null'
      - name: interval_merging_rule
        symbols: [ALL, OVERLAPPING_ONLY]
        type: enum
    - doc: Amount of padding (in bp) to add to each interval you are including.
      id: interval_padding
      inputBinding: {position: 4, prefix: --interval-padding, shellQuote: false}
      label: Interval padding
      sbg:altPrefix: -ip
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: Set merging approach to use for combining interval inputs.
      id: interval_set_rule
      inputBinding: {position: 4, prefix: --interval-set-rule, shellQuote: false}
      label: Interval set rule
      sbg:altPrefix: -isr
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: UNION
      type:
      - 'null'
      - name: interval_set_rule
        symbols: [UNION, INTERSECTION]
        type: enum
    - doc: File which contains one or more genomic intervals over which to operate.
      id: include_intervals_file
      inputBinding: {position: 4, prefix: --intervals, shellQuote: false}
      label: Include intervals file
      sbg:altPrefix: -L
      sbg:category: Optional Arguments
      sbg:fileTypes: BED
      sbg:toolDefaultValue: 'null'
      type: ['null', File]
    - doc: One or more genomic intervals over which to operate.
      id: include_intervals_string
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--intervals', self[i]);\n
          \       }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Include intervals string
      sbg:altPrefix: -L
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: The name of the read group to keep. Valid only if "ReadGroupReadFilter"
        is specified.
      id: keep_read_group
      inputBinding: {position: 4, prefix: --keep-read-group, shellQuote: false}
      label: Keep read group
      sbg:category: Conditional Arguments
      type: ['null', string]
    - doc: Keep only reads on the reverse strand. Valid only if "ReadStrandFilter"
        is specified
      id: keep_reverse_strand_only
      inputBinding: {position: 4, prefix: --keep-reverse-strand-only, shellQuote: false}
      label: Keep reverse strand only
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Name of the library to keep. Valid only if "LibraryReadFilter" is specified
      id: library
      inputBinding: {position: 4, shellQuote: false, valueFrom: "${\n    if (self)\n
          \   {\n        var cmd = [];\n        for (var i = 0; i < self.length; i++)
          \n        {\n            cmd.push('--library', self[i]);\n        }\n        return
          cmd.join(' ');\n    }\n    \n}"}
      label: Library
      sbg:altPrefix: -library
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Maximum length of fragment (insert size). Valid only if "FragmentLengthReadFilter"
        is specified.
      id: max_fragment_length
      inputBinding: {position: 4, prefix: --max-fragment-length, shellQuote: false}
      label: Max fragment length
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '1000000'
      type: ['null', int]
    - doc: Keep only reads with length at most equal to the specified value. Valid
        only if "ReadLengthReadFilter" is specified.
      id: max_read_length
      inputBinding: {position: 4, prefix: --max-read-length, shellQuote: false}
      label: Max read length
      sbg:category: Conditional Arguments
      type: ['null', int]
    - doc: Maximum mapping quality to keep (inclusive). Valid only if "MappingQualityReadFilter"
        is specified.
      id: maximum_mapping_quality
      inputBinding: {position: 4, prefix: --maximum-mapping-quality, shellQuote: false}
      label: Maximum mapping quality
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', int]
    - doc: This input allows a user to set the desired overhead memory when running
        a tool or adding it to a workflow. This amount will be added to the Memory
        per job in the Memory requirements section but it will not be added to the
        -Xmx parameter leaving some memory not occupied which can be used as stack
        memory (-Xmx parameter defines heap memory). This input should be defined
        in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      id: memory_overhead_per_job
      label: Memory overhead per job
      sbg:category: Platform Options
      type: ['null', int]
    - doc: This input allows a user to set the desired memory requirement when running
        a tool or adding it to a workflow. This value should be propagated to the
        -Xmx parameter too.This input should be defined in MB (for both the platform
        part and the -Xmx part if Java tool is wrapped).
      id: memory_per_job
      label: Memory per job
      sbg:category: Platform Options
      type: ['null', int]
    - doc: Keep only reads with length at least equal to the specified value. Valid
        only if "ReadLengthReadFilter" is specified.
      id: min_read_length
      inputBinding: {position: 4, prefix: --min-read-length, shellQuote: false}
      label: Min read length
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    - doc: Minimum mapping quality to keep (inclusive). Valid only if "MappingQualityReadFilter"
        is specified.
      id: minimum_mapping_quality
      inputBinding: {position: 4, prefix: --minimum-mapping-quality, shellQuote: false}
      label: Minimum mapping quality
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '10'
      type: ['null', int]
    - doc: Platform attribute (PL) to match. Valid only if "PlatformReadFilter" is
        specified.
      id: platform_filter_name
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--platform-filter-name',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Platform filter name
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Don't recalibrate bases with quality scores less than this threshold.
      id: preserve_qscores_less_than
      inputBinding: {position: 4, prefix: --preserve-qscores-less-than, shellQuote: false}
      label: Preserve qscores less than
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '6'
      type: ['null', int]
    - doc: Quantize quality scores to a given number of levels. Cannot be used in
        conjuction with argument(s) staticQuantizationQuals, roundDown.
      id: quantize_quals
      inputBinding: {position: 4, prefix: --quantize-quals, shellQuote: false}
      label: Quantize quals
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: Read filters to be applied before analysis.
      id: read_filter
      inputBinding: {position: 4, shellQuote: false, valueFrom: "${\n    if (self)\n
          \   {\n        var cmd = [];\n        for (var i = 0; i < self.length; i++)
          \n        {\n            cmd.push('--read-filter', self[i]);\n        }\n
          \       return cmd.join(' ');\n    }\n    \n}"}
      label: Read filter
      sbg:altPrefix: -RF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - items:
          name: read_filter
          symbols: [AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter,
            AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter,
            FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter,
            LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter,
            MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter,
            MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter,
            MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter,
            NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter,
            NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter,
            OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter,
            PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter,
            ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter,
            ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter,
            SecondOfPairReadFilter, SeqIsStoredReadFilter, ValidAlignmentEndReadFilter,
            ValidAlignmentStartReadFilter, WellformedReadFilter]
          type: enum
        type: array
    - doc: The name of the read group to filter out. Valid only if "ReadGroupBlackListReadFilter"
        is specified.
      id: read_group_black_list
      inputBinding: {position: 4, shellQuote: false, valueFrom: "${\n    if (self)\n
          \   {\n        var cmd = [];\n        for (var i = 0; i < self.length; i++)
          \n        {\n            cmd.push('--read-group-black-list', self[i]);\n
          \       }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Read group black list
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Keep only reads with this read name. Valid only if "ReadNameReadFilter"
        is specified.
      id: read_name
      inputBinding: {position: 4, prefix: --read-name, shellQuote: false}
      label: Read name
      sbg:category: Conditional Arguments
      type: ['null', string]
    - doc: Validation stringency for all SAM/BAM/CRAM files read by this program.
        The default stringency value SILENT can improve performance when processing
        a BAM file in which variable-length data (read, qualities, tags) do not otherwise
        need to be decoded.
      id: read_validation_stringency
      inputBinding: {position: 4, prefix: --read-validation-stringency, shellQuote: false}
      label: Read validation stringency
      sbg:altPrefix: -VS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: SILENT
      type:
      - 'null'
      - name: read_validation_stringency
        symbols: [STRICT, LENIENT, SILENT]
        type: enum
    - doc: Reference sequence.
      id: in_reference
      inputBinding: {position: 4, prefix: --reference, shellQuote: false}
      label: Reference
      sbg:altPrefix: -R
      sbg:category: Optional Arguments
      sbg:fileTypes: FASTA, FA
      sbg:toolDefaultValue: 'null'
      secondaryFiles: [.fai, ^.dict]
      type: ['null', File]
    - doc: Round quals down to nearest quantized qual. Cannot be used in conjuction
        with argument quantizationLevels.
      id: round_down_quantized
      inputBinding: {position: 4, prefix: --round-down-quantized, shellQuote: false}
      label: Round down quantized
      sbg:category: Advanced Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: The name of the sample(s) to keep, filtering out all others. Valid only
        if "SampleReadFilter" is specified.
      id: sample
      inputBinding: {position: 4, shellQuote: false, valueFrom: "${\n    if (self)\n
          \   {\n        var cmd = [];\n        for (var i = 0; i < self.length; i++)
          \n        {\n            cmd.push('--sample', self[i]);\n        }\n        return
          cmd.join(' ');\n    }\n    \n}"}
      label: Sample
      sbg:altPrefix: -sample
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Output traversal statistics every time this many seconds elapse.
      id: seconds_between_progress_updates
      inputBinding: {position: 4, prefix: --seconds-between-progress-updates, shellQuote: false}
      label: Seconds between progress updates
      sbg:altPrefix: -seconds-between-progress-updates
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '10.0'
      type: ['null', float]
    - doc: Use the given sequence dictionary as the master/canonical sequence dictionary.
        Must be a .dict file.
      id: sequence_dictionary
      inputBinding: {position: 4, prefix: --sequence-dictionary, shellQuote: false}
      label: Sequence dictionary
      sbg:altPrefix: -sequence-dictionary
      sbg:category: Optional Arguments
      sbg:fileTypes: DICT
      sbg:toolDefaultValue: '10.0'
      type: ['null', File]
    - doc: Use static quantized quality scores to a given number of levels (with -bqsr).
        Cannot be used in conjuction with argument(s) quantizationLevels.
      id: static_quantized_quals
      inputBinding: {position: 4, shellQuote: false, valueFrom: "${\n    if (self)\n
          \   {\n        var cmd = [];\n        for (var i = 0; i < self.length; i++)
          \n        {\n            cmd.push('--static-quantized-quals', self[i]);\n
          \       }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Static quantized quals
      sbg:category: Advanced Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - {items: int, type: array}
    - doc: Use the base quality scores from the OQ tag.
      id: use_original_qualities
      inputBinding: {position: 4, prefix: --use-original-qualities, shellQuote: false}
      label: Use original qualities
      sbg:altPrefix: -OQ
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: This input allows a user to set the desired CPU requirement when running
        a tool or adding it to a workflow.
      id: cpu_per_job
      label: CPU per job
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    - doc: Output file name prefix.
      id: output_prefix
      label: Output prefix
      sbg:category: Optional Arguments
      type: ['null', string]
    - doc: Output file format
      id: output_file_format
      label: Output file format
      sbg:category: Optional Arguments
      type:
      - 'null'
      - name: output_file_format
        symbols: [sam, bam, cram]
        type: enum
    label: GATK ApplyBQSR
    outputs:
    - doc: Output base quality recalibrated BAM/SAM/CRAM file.
      id: out_alignments
      label: Output recalibrated BAM/SAM/CRAM file
      outputBinding: {glob: '*am', outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: BAM, SAM, CRAM
      secondaryFiles: ["${\n    var output_ext = self.nameext;\n    \n    if (output_ext
          == \".bam\")\n    {\n        return self.nameroot + \".bai\";\n    }\n    else
          if (output_ext == \".cram\")\n    {\n        return self.nameroot + \".crai\";\n
          \   }\n    else \n        return ''; \n  \n}"]
      type: ['null', File]
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: "${\n    return inputs.cpu_per_job ?
        inputs.cpu_per_job : 1;\n}", ramMin: "${\n    var memory = 3500;\n    if (inputs.memory_per_job)
        \n    {\n        memory = inputs.memory_per_job;\n    }\n    if (inputs.memory_overhead_per_job)\n
        \   {\n        memory += inputs.memory_overhead_per_job;\n    }\n    return
        memory;\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n
          \   else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar
          groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict
          = {};\n    for (var i = 0; i < files.length; i++) {\n        var value =
          files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};", "\nvar setMetadata = function(file,
          metadata) {\n    if (!('metadata' in file))\n        file['metadata'] =
          metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BAM Processing]
    sbg:content_hash: a0fe80d5ce988982acb05bdf82d9c9874f2ec68e4160fc40b157059decbfb5810
    sbg:contributors: [uros_sipetic, nemanja.vucic, nens, veliborka_josipovic]
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/24
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552923344
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-applybqsr-4-1-0-0/10
    sbg:image_url: null
    sbg:latestRevision: 10
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - {id: 'https://software.broadinstitute.org/gatk/', label: Homepage}
    - {id: 'https://github.com/broadinstitute/gatk/', label: Source Code}
    - {id: 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip',
      label: Download}
    - {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
    - {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php',
      label: Documentation}
    sbg:modifiedBy: nemanja.vucic
    sbg:modifiedOn: 1559744828
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 10
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/24
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1552923344, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/8}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554493022, 'sbg:revision': 1,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/14}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554493059, 'sbg:revision': 2,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/15}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554720859, 'sbg:revision': 3,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/16}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554999197, 'sbg:revision': 4,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/17}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557734544, 'sbg:revision': 5, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/18}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000590, 'sbg:revision': 6, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/19}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351541, 'sbg:revision': 7, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/21}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558451164, 'sbg:revision': 8, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/22}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558524331, 'sbg:revision': 9, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/23}
    - {'sbg:modifiedBy': nemanja.vucic, 'sbg:modifiedOn': 1559744828, 'sbg:revision': 10,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-applybqsr-4-1-0-0/24}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: 1615.560546875
  sbg:y: 207.82618713378906
  scatter: [include_intervals_file]
- id: gatk_baserecalibrator_4_1_0_0
  in:
  - {id: in_alignments, source: gatk_setnmmdanduqtags_4_1_0_0/out_alignments}
  - {id: include_intervals_file, source: sbg_lines_to_interval_list/out_intervals}
  - id: known_indels
    source: [known_indels]
  - id: known_snps
    source: [known_snps]
  - {id: in_reference, source: in_reference}
  - {default: true, id: use_original_qualities}
  label: GATK BaseRecalibrator
  out:
  - {id: output}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, prefix: '', shellQuote: false, valueFrom: /opt/gatk}
    - {position: 1, shellQuote: false, valueFrom: --java-options}
    - {position: 2, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)
        {\n        return '\\\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\\\"';\n
        \   }\n    return '\\\"-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal
        -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log
        -Xms4000m\\\"';\n}"}
    - {position: 3, shellQuote: false, valueFrom: BaseRecalibrator}
    - {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n    var in_alignments
        = inputs.in_alignments;\n    var output_prefix = '';\n    if (inputs.output_prefix)\n
        \   {\n        output_prefix = inputs.output_prefix;\n    }\n    else \n    {\n
        \       if (in_alignments.metadata && in_alignments.metadata.sample_id)\n
        \       {\n             output_prefix = in_alignments.metadata.sample_id;\n
        \       }\n        else \n        {\n            output_prefix = in_alignments.path.split('/').pop().split('.')[0];\n
        \       }\n    }\n    \n    return \"--output \" + output_prefix + \".recal_data.csv\";\n}"}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "The **GATK BaseRecalibrator** tool performs first pass of Base Quality Score
      Recalibration (BQSR) based on various covariates. The default covariates are
      read group, reported quality score, machine cycle, and nucleotide context. The
      tool generates on a recalibration table on its output.\n\n*A list of **all inputs
      and parameters** with corresponding descriptions can be found at the bottom
      of the page.*\n\n###Common Use Cases\n\n* The **GATK BaseRecalibrator** tool
      requires the input read data whose quality scores need to be assessed on its
      **Input BAM/SAM/CRAM file** (`--input`) input and the database of known polymorphic
      sites to skip over on its **Known SNPs** and **Known INDELs** (`--known-sites`)
      inputs. The tool generates on its **Output recalibration report** output a GATK
      report file with many tables: the list of arguments, the quantized qualities
      table, the recalibration table by read group, the recalibration table by quality
      score,\nthe recalibration table for all the optional covariates.\n\n* Usage
      example:\n\n```\ngatk BaseRecalibrator \\\n   --input my_reads.bam \\\n   --reference
      reference.fasta \\\n   --known-sites sites_of_variation.vcf \\\n   --known-sites
      another/optional/setOfSitesToMask.vcf \\\n   --output recal_data.table\n\n```\n\n###Changes
      Introduced by Seven Bridges\n\n* All output files will be prefixed using the
      **Output prefix** parameter. In case **Output prefix** is not provided, output
      prefix will be the same as the Sample ID metadata from **Input SAM/BAM/CRAM
      file**, if the Sample ID metadata exists. Otherwise, output prefix will be inferred
      from the **Input SAM/BAM/CRAM** filename. This way, having identical names of
      the output files between runs is avoided. Moreover,  **recal_data** will be
      added before the extension of the output file name which is CSV by default.
      \n\n* **Include intervals** (`--intervals`) option is divided into **Include
      intervals string** and **Include intervals file** options.\n\n* **Exclude intervals**
      (`--exclude-intervals`) option is divided into **Exclude intervals string**
      and **Exclude intervals file** options.\n\n###Common Issues and Important Notes\n\n*
      Note: If **Read filter** (`--read-filter`) option is set to \"LibraryReadFilter\",
      **Library** (`--library`) option must be set to some value.\n* Note: If **Read
      filter** (`--read-filter`) option is set to \"PlatformReadFilter\", **Platform
      filter name** (`--platform-filter-name`) option must be set to some value.\n*
      Note: If **Read filter** (`--read-filter`) option is set to\"PlatformUnitReadFilter\",
      **Black listed lanes** (`--black-listed-lanes`) option must be set to some value.
      \n* Note: If **Read filter** (`--read-filter`) option is set to \"ReadGroupBlackListReadFilter\",
      **Read group black list** (`--read-group-black-list`) option must be set to
      some value.\n* Note: If **Read filter** (`--read-filter`) option is set to \"ReadGroupReadFilter\",
      **Keep read group** (`--keep-read-group`) option must be set to some value.\n*
      Note: If **Read filter** (`--read-filter`) option is set to \"ReadLengthReadFilter\",
      **Max read length** (`--max-read-length`) option must be set to some value.\n*
      Note: If **Read filter** (`--read-filter`) option is set to \"ReadNameReadFilter\",
      **Read name** (`--read-name`) option must be set to some value.\n* Note: If
      **Read filter** (`--read-filter`) option is set to \"ReadStrandFilter\", **Keep
      reverse strand only** (`--keep-reverse-strand-only`) option must be set to some
      value.\n* Note: If **Read filter** (`--read-filter`) option is set to \"SampleReadFilter\",
      **Sample** (`--sample`) option must be set to some value.\n\n###Performance
      Benchmarking\n\nBelow is a table describing runtimes and task costs of **GATK
      BaseRecalibrator** for a couple of different samples, executed on the AWS cloud
      instances:\n\n| Experiment type |  Input size | Duration |  Cost | Instance
      (AWS) | \n|:--------------:|:------------:|:--------:|:-------:|:---------:|\n|
      \    RNA-Seq     |  2.2 GB |   8min   | ~0.05$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq
      \    |  6.6 GB |   18min   | ~0.12$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq
      \    | 11 GB |  26min  | ~0.17$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq     |
      22 GB |  45min  | ~0.29$ | c4.2xlarge (8 CPUs) |\n\n*Cost can be significantly
      reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances)
      for more details.*\n\n###References\n\n[1] [GATK BaseRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php)"
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-baserecalibrator-4-1-0-0/11
    inputs:
    - doc: Threshold number of ambiguous bases. If null, uses threshold fraction;
        otherwise, overrides threshold fraction. Cannot be used in conjuction with
        argument(s) maxAmbiguousBaseFraction. Valid only if "AmbiguousBaseReadFilter"
        is specified.
      id: ambig_filter_bases
      inputBinding: {position: 4, prefix: --ambig-filter-bases, shellQuote: false}
      label: Ambig filter bases
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', int]
    - doc: Threshold fraction of ambiguous bases. Cannot be used in conjuction with
        argument(s) maxAmbiguousBases. Valid only if "AmbiguousBaseReadFilter" is
        specified.
      id: ambig_filter_frac
      inputBinding: {position: 4, prefix: --ambig-filter-frac, shellQuote: false}
      label: Ambig filter frac
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', float]
    - doc: The binary tag covariate name if using it.
      id: binary_tag_name
      inputBinding: {position: 4, prefix: --binary-tag-name, shellQuote: false}
      label: Binary tag name
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: Platform unit (PU) to filter out. This argument must be specified at least
        once. Valid only if "PlatformUnitReadFilter" is specified.
      id: black_listed_lanes
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--black-listed-lanes',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Black listed lanes
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: BQSR BAQ gap open penalty (Phred Scaled). Default value is 40. 30 is perhaps
        better for whole genome call sets.
      id: bqsr_baq_gap_open_penalty
      inputBinding: {position: 4, prefix: --bqsr-baq-gap-open-penalty, shellQuote: false}
      label: BQSR BAQ gap open penalty
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '40'
      type: ['null', float]
    - doc: Assign a default base quality.
      id: default_base_qualities
      inputBinding: {position: 4, prefix: --default-base-qualities, shellQuote: false}
      label: Default base qualities
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '-1'
      type: ['null', int]
    - doc: Default quality for the base deletions covariate.
      id: deletions_default_quality
      inputBinding: {position: 4, prefix: --deletions-default-quality, shellQuote: false}
      label: Deletions default quality
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '45'
      type: ['null', int]
    - doc: Read filters to be disabled before analysis.
      id: disable_read_filter
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--disable-read-filter',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Disable read filter
      sbg:altPrefix: -DF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - items:
          name: disable_read_filter
          symbols: [AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter,
            AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter,
            FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter,
            LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter,
            MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter,
            MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter,
            MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter,
            NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter,
            NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter,
            OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter,
            PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter,
            ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter,
            ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter,
            SecondOfPairReadFilter, SeqIsStoredReadFilter, ValidAlignmentEndReadFilter,
            ValidAlignmentStartReadFilter, WellformedReadFilter]
          type: enum
        type: array
    - doc: If specified, do not check the sequence dictionaries from our inputs for
        compatibility. Use at your own risk!
      id: disable_sequence_dictionary_validation
      inputBinding: {position: 4, prefix: --disable-sequence-dictionary-validation,
        shellQuote: false}
      label: Disable sequence dictionary validation
      sbg:altPrefix: -disable-sequence-dictionary-validation
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: 'Disable all tool default read filters (WARNING: many tools will not function
        correctly without their default read filters on).'
      id: disable_tool_default_read_filters
      inputBinding: {position: 4, prefix: --disable-tool-default-read-filters, shellQuote: false}
      label: Disable tool default read filters
      sbg:altPrefix: -disable-tool-default-read-filters
      sbg:category: Advanced Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Allow a read to be filtered out based on having only 1 soft-clipped block.
        By default, both ends must have a soft-clipped block, setting this flag requires
        only 1 soft-clipped block. Valid only if "OverclippedReadFilter" is specified.
      id: dont_require_soft_clips_both_ends
      inputBinding: {position: 4, prefix: --dont-require-soft-clips-both-ends, shellQuote: false}
      label: Dont require soft clips both ends
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: File which contains one or more genomic intervals to exclude from processing.
      id: exclude_intervals_file
      inputBinding: {position: 4, prefix: --exclude-intervals, shellQuote: false}
      label: Exclude intervals file
      sbg:altPrefix: -XL
      sbg:category: Optional Arguments
      sbg:fileTypes: BED
      sbg:toolDefaultValue: 'null'
      type: ['null', File]
    - doc: One or more genomic intervals to exclude from processing.
      id: exclude_intervals_string
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--exclude-intervals',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Exclude intervals string
      sbg:altPrefix: -XL
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Minimum number of aligned bases. Valid only if "OverclippedReadFilter"
        is specified.
      id: filter_too_short
      inputBinding: {position: 4, prefix: --filter-too-short, shellQuote: false}
      label: Filter too short
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '30'
      type: ['null', int]
    - doc: Size of the k-mer context to be used for base insertions and deletions.
      id: indels_context_size
      inputBinding: {position: 4, prefix: --indels-context-size, shellQuote: false}
      label: Indels context size
      sbg:altPrefix: -ics
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '3'
      type: ['null', int]
    - doc: BAM/SAM/CRAM file containing reads.
      id: in_alignments
      inputBinding: {position: 4, prefix: --input, shellQuote: false}
      label: Input BAM/SAM/CRAM file
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM, CRAM
      secondaryFiles: ["${\n    if (self.nameext == \".bam\")\n    {\n        return
          \ self.nameroot + \".bai\";\n    }\n    else if (self.nameext == \".cram\")\n
          \   {\n        return  self.nameroot + \".crai\";\n    }\n    return '';
          \n    \n}"]
      type: File
    - doc: Default quality for the base insertions covariate.
      id: insertions_default_quality
      inputBinding: {position: 4, prefix: --insertions-default-quality, shellQuote: false}
      label: Insertions default quality
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '45'
      type: ['null', int]
    - doc: Amount of padding (in bp) to add to each interval you are excluding.
      id: interval_exclusion_padding
      inputBinding: {position: 4, prefix: --interval-exclusion-padding, shellQuote: false}
      label: Interval exclusion padding
      sbg:altPrefix: -ixp
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: Interval merging rule for abutting intervals.
      id: interval_merging_rule
      inputBinding: {position: 4, prefix: --interval-merging-rule, shellQuote: false}
      label: Interval merging rule
      sbg:altPrefix: -imr
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: ALL
      type:
      - 'null'
      - name: interval_merging_rule
        symbols: [ALL, OVERLAPPING_ONLY]
        type: enum
    - doc: Amount of padding (in bp) to add to each interval you are including.
      id: interval_padding
      inputBinding: {position: 4, prefix: --interval-padding, shellQuote: false}
      label: Interval padding
      sbg:altPrefix: -ip
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: Set merging approach to use for combining interval inputs.
      id: interval_set_rule
      inputBinding: {position: 4, prefix: --interval-set-rule, shellQuote: false}
      label: Interval set rule
      sbg:altPrefix: -isr
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: UNION
      type:
      - 'null'
      - name: interval_set_rule
        symbols: [UNION, INTERSECTION]
        type: enum
    - doc: File which contains one or more genomic intervals over which to operate.
      id: include_intervals_file
      inputBinding: {position: 4, prefix: --intervals, shellQuote: false}
      label: Include intervals file
      sbg:altPrefix: -L
      sbg:category: Optional Arguments
      sbg:fileTypes: BED, INTERVALS
      sbg:toolDefaultValue: 'null'
      type: ['null', File]
    - doc: One or more genomic intervals over which to operate.
      id: include_intervals_string
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--intervals', self[i]);\n
          \       }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Include intervals string
      sbg:altPrefix: -L
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: The name of the read group to keep. Valid only if "ReadGroupReadFilter"
        is specified.
      id: keep_read_group
      inputBinding: {position: 4, prefix: --keep-read-group, shellQuote: false}
      label: Keep read group
      sbg:category: Conditional Arguments
      type: ['null', string]
    - doc: Keep only reads on the reverse strand. Valid only if "ReadStrandFilter"
        is specified.
      id: keep_reverse_strand_only
      inputBinding: {position: 4, prefix: --keep-reverse-strand-only, shellQuote: false}
      label: Keep reverse strand only
      sbg:category: Conditional Arguments
      type: ['null', boolean]
    - doc: One or more databases of known polymorphic sites used to exclude regions
        around known polymorphisms from analysis. This argument must be specified
        at least once.
      id: known_indels
      inputBinding: {position: 5, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--known-sites', self[i].path);\n
          \       }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Known indels
      sbg:category: Required Arguments
      sbg:fileTypes: VCF
      secondaryFiles: ["${\n    return self.basename + \".idx\"\n}"]
      type:
      - 'null'
      - {items: File, type: array}
    - doc: One or more databases of known polymorphic sites used to exclude regions
        around known polymorphisms from analysis. This argument must be specified
        at least once.
      id: known_snps
      inputBinding: {position: 5, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--known-sites', self[i].path);\n
          \       }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Known SNPs
      sbg:category: Required Arguments
      sbg:fileTypes: VCF
      secondaryFiles: ["${\n    return self.basename + \".idx\"\n}"]
      type:
      - 'null'
      - {items: File, type: array}
    - doc: Name of the library to keep. This argument must be specified at least once.
        Valid only if "LibraryReadFilter" is specified.
      id: library
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--library', self[i]);\n
          \       }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Library
      sbg:altPrefix: -library
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Minimum quality for the bases in the tail of the reads to be considered.
      id: low_quality_tail
      inputBinding: {position: 4, prefix: --low-quality-tail, shellQuote: false}
      label: Low quality tail
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
      type: ['null', int]
    - doc: Maximum length of fragment (insert size). Valid only if "FragmentLengthReadFilter"
        is specified.
      id: max_fragment_length
      inputBinding: {position: 4, prefix: --max-fragment-length, shellQuote: false}
      label: Max fragment length
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '1000000'
      type: ['null', int]
    - doc: Keep only reads with length at most equal to the specified value. Valid
        only if "ReadLengthReadFilter" is specified.
      id: max_read_length
      inputBinding: {position: 4, prefix: --max-read-length, shellQuote: false}
      label: Max read length
      sbg:category: Conditional Arguments
      type: ['null', int]
    - doc: The maximum cycle value permitted for the Cycle covariate.
      id: maximum_cycle_value
      inputBinding: {position: 4, prefix: --maximum-cycle-value, shellQuote: false}
      label: Maximum cycle value
      sbg:altPrefix: -max-cycle
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500'
      type: ['null', int]
    - doc: Maximum mapping quality to keep (inclusive). Valid only if "MappingQualityReadFilter"
        is specified.
      id: maximum_mapping_quality
      inputBinding: {position: 4, prefix: --maximum-mapping-quality, shellQuote: false}
      label: Maximum mapping quality
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', int]
    - doc: This input allows a user to set the desired overhead memory when running
        a tool or adding it to a workflow. This amount will be added to the Memory
        per job in the Memory requirements section but it will not be added to the
        -Xmx parameter leaving some memory not occupied which can be used as stack
        memory (-Xmx parameter defines heap memory). This input should be defined
        in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      id: memory_overhead_per_job
      label: Memory overhead per job
      sbg:category: Platform Options
      type: ['null', int]
    - doc: This input allows a user to set the desired memory requirement when running
        a tool or adding it to a workflow. This value should be propagated to the
        -Xmx parameter too.This input should be defined in MB (for both the platform
        part and the -Xmx part if Java tool is wrapped).
      id: memory_per_job
      label: Memory per job
      sbg:category: Platform Options
      type: ['null', int]
    - doc: Keep only reads with length at least equal to the specified value. Valid
        only if "ReadLengthReadFilter" is specified.
      id: min_read_length
      inputBinding: {position: 4, prefix: --min-read-length, shellQuote: false}
      label: Min read length
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    - doc: Minimum mapping quality to keep (inclusive). Valid only if "MappingQualityReadFilter"
        is specified.
      id: minimum_mapping_quality
      inputBinding: {position: 4, prefix: --minimum-mapping-quality, shellQuote: false}
      label: Minimum mapping quality
      sbg:category: Conditional Arguments
      sbg:toolDefaultValue: '10'
      type: ['null', int]
    - doc: Size of the k-mer context to be used for base mismatches.
      id: mismatches_context_size
      inputBinding: {position: 4, prefix: --mismatches-context-size, shellQuote: false}
      label: Mismatches context size
      sbg:altPrefix: -mcs
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
      type: ['null', int]
    - doc: Default quality for the base mismatches covariate.
      id: mismatches_default_quality
      inputBinding: {position: 4, prefix: --mismatches-default-quality, shellQuote: false}
      label: Mismatches default quality
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '-1'
      type: ['null', int]
    - doc: Platform attribute (PL) to match. This argument must be specified at least
        once. Valid only if "PlatformReadFilter" is specified.
      id: platform_filter_name
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--platform-filter-name',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Platform filter name
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Don't recalibrate bases with quality scores less than this threshold (with
        -bqsr).
      id: preserve_qscores_less_than
      inputBinding: {position: 4, prefix: --preserve-qscores-less-than, shellQuote: false}
      label: Preserve qscores less than
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '6'
      type: ['null', int]
    - doc: Number of distinct quality scores in the quantized output.
      id: quantizing_levels
      inputBinding: {position: 4, prefix: --quantizing-levels, shellQuote: false}
      label: Quantizing levels
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '16'
      type: ['null', int]
    - doc: Read filters to be applied before analysis.
      id: read_filter
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--read-filter', self[i]);\n
          \       }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Read filter
      sbg:altPrefix: -RF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - items:
          name: read_filter
          symbols: [AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter,
            AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter,
            FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter,
            LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter,
            MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter,
            MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter,
            MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter,
            NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter,
            NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter,
            OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter,
            PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter,
            ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter,
            ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter,
            SecondOfPairReadFilter, SeqIsStoredReadFilter, ValidAlignmentEndReadFilter,
            ValidAlignmentStartReadFilter, WellformedReadFilter]
          type: enum
        type: array
    - doc: The name of the read group to filter out. This argument must be specified
        at least once. Valid only if "ReadGroupBlackListReadFilter" is specified.
      id: read_group_black_list
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--read-group-black-list',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Read group black list
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Keep only reads with this read name. Valid only if "ReadNameReadFilter"
        is specified.
      id: read_name
      inputBinding: {position: 4, prefix: --read-name, shellQuote: false}
      label: Read name
      sbg:category: Conditional Arguments
      type: ['null', string]
    - doc: Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.
        The default stringency value SILENT can improve performance when processing
        a BAM file in which variable-length data (read, qualities, tags) do not otherwise
        need to be decoded.
      id: read_validation_stringency
      inputBinding: {position: 4, prefix: --read-validation-stringency, shellQuote: false}
      label: Read validation stringency
      sbg:altPrefix: -VS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: SILENT
      type:
      - 'null'
      - name: read_validation_stringency
        symbols: [STRICT, LENIENT, SILENT]
        type: enum
    - doc: Reference sequence file.
      id: in_reference
      inputBinding: {position: 4, prefix: --reference, shellQuote: false}
      label: Reference
      sbg:altPrefix: -R
      sbg:category: Required Arguments
      sbg:fileTypes: FASTA, FA
      secondaryFiles: [.fai, ^.dict]
      type: File
    - doc: The name of the sample(s) to keep, filtering out all others. This argument
        must be specified at least once. Valid only if "SampleReadFilter" is specified.
      id: sample
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--sample', self[i]);\n
          \       }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Sample
      sbg:altPrefix: -sample
      sbg:category: Conditional Arguments
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Use the given sequence dictionary as the master/canonical sequence dictionary.
        Must be a .dict file.
      id: sequence_dictionary
      inputBinding: {position: 4, prefix: --sequence-dictionary, shellQuote: false}
      label: Sequence dictionary
      sbg:altPrefix: -sequence-dictionary
      sbg:category: Optional Arguments
      sbg:fileTypes: DICT
      sbg:toolDefaultValue: '10.0'
      type: ['null', File]
    - doc: Use the base quality scores from the OQ tag.
      id: use_original_qualities
      inputBinding: {position: 4, prefix: --use-original-qualities, shellQuote: false}
      label: Use original qualities
      sbg:altPrefix: -OQ
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Output file name prefix.
      id: output_prefix
      label: Output prefix
      sbg:category: Optional arguments
      type: ['null', string]
    - doc: This input allows a user to set the desired CPU requirement when running
        a tool or adding it to a workflow.
      id: cpu_per_job
      label: CPU per job
      sbg:category: Platform Options
      type: ['null', int]
    label: GATK BaseRecalibrator
    outputs:
    - doc: The output recalibration table file.
      id: output
      label: Output recalibration report
      outputBinding: {glob: '*csv', outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: CSV
      type: ['null', File]
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: "${\n    return inputs.cpu_per_job ?
        inputs.cpu_per_job : 1;\n}", ramMin: "${\n    var memory = 6500;\n    if (inputs.memory_per_job)
        \n    {\n        memory = inputs.memory_per_job;\n    }\n    if (inputs.memory_overhead_per_job)\n
        \   {\n        memory += inputs.memory_overhead_per_job;\n    }\n    return
        memory;\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n
          \   else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar
          groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict
          = {};\n    for (var i = 0; i < files.length; i++) {\n        var value =
          files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};", "\nvar setMetadata = function(file,
          metadata) {\n    if (!('metadata' in file))\n        file['metadata'] =
          metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BAM Processing]
    sbg:content_hash: af9bd2f52de83a0244fc631865ff864f27bcc1a64126b96ea1b89e70f2e78cd4a
    sbg:contributors: [nens, uros_sipetic, veliborka_josipovic]
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/26
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552922094
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-baserecalibrator-4-1-0-0/11
    sbg:image_url: null
    sbg:latestRevision: 11
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - {id: 'https://software.broadinstitute.org/gatk/', label: Homepage}
    - {id: 'https://github.com/broadinstitute/gatk/', label: Source Code}
    - {id: 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip',
      label: Download}
    - {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
    - {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php',
      label: Documentation}
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1558518057
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 11
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/26
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1552922094, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/11}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554492924, 'sbg:revision': 1,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/14}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554492998, 'sbg:revision': 2,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/15}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554720866, 'sbg:revision': 3,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/17}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554999207, 'sbg:revision': 4,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/18}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1556030757, 'sbg:revision': 5,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/19}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557735256, 'sbg:revision': 6, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/20}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000594, 'sbg:revision': 7, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/21}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351546, 'sbg:revision': 8, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/23}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558450805, 'sbg:revision': 9, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/24}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558517350, 'sbg:revision': 10, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/25}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558518057, 'sbg:revision': 11, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-baserecalibrator-4-1-0-0/26}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: 1241.2686767578125
  sbg:y: 307.5648193359375
  scatter: [include_intervals_file]
- id: gatk_createsequencegroupingtsv_4_1_0_0
  in:
  - {id: ref_dict, source: ref_dict}
  label: GATK CreateSequenceGroupingTSV
  out:
  - {id: sequence_grouping}
  - {id: sequence_grouping_with_unmapped}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    baseCommand: [python, CreateSequenceGroupingTSV.py]
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "**CreateSequenceGroupingTSV** tool generate sets of intervals for scatter-gathering
      over chromosomes.\n\nIt takes **Reference dictionary** file (`--ref_dict`) as
      an input and creates files which contain chromosome names grouped based on their
      sizes.\n\n\n###**Common Use Cases**\n\nThe tool has only one input (`--ref_dict`)
      which is required and has no additional arguments. **CreateSequenceGroupingTSV**
      tool results are **Sequence Grouping** file which is a text file containing
      chromosome groups, and **Sequence Grouping with Unmapped**, a text file which
      has the same content as **Sequence Grouping** with additional line containing
      \"unmapped\" string.\n\n\n* Usage example\n\n\n```\npython CreateSequenceGroupingTSV.py
      \n      --ref_dict example_reference.dict\n\n```\n\n\n\n###**Changes Introduced
      by Seven Bridges**\n\nPython code provided within WGS Germline WDL was adjusted
      to be called as a script (`CreateSequenceGroupingTSV.py`).\n\n\n###**Common
      Issues and Important Notes**\n\nNone.\n\n\n### Reference\n[1] [CreateSequenceGroupingTSV](https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf-fc-hg38.wdl)"
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-createsequencegroupingtsv-4-1-0-0/4
    inputs:
    - doc: Reference dictionary containing information about chromosome names and
        their lengths.
      id: ref_dict
      inputBinding: {position: 0, prefix: --ref_dict, shellQuote: false}
      label: Reference Dictionary
      sbg:fileTypes: DICT
      type: File
    label: GATK CreateSequenceGroupingTSV
    outputs:
    - doc: Each line of the file represents one group of chromosomes which are processed
        together in later steps of the GATK Germline workflow. The groups are determined
        based on the chromosomes sizes.
      id: sequence_grouping
      label: Sequence Grouping
      outputBinding: {glob: sequence_grouping.txt, outputEval: '$(inheritMetadata(self,
          inputs.ref_dict))'}
      sbg:fileTypes: TXT
      type: ['null', File]
    - doc: The file has the same content as "Sequence Grouping" file, with an additional,
        last line containing "unmapped" string.
      id: sequence_grouping_with_unmapped
      label: Sequence Grouping with Unmapped
      outputBinding: {glob: sequence_grouping_with_unmapped.txt, outputEval: '$(inheritMetadata(self,
          inputs.ref_dict))'}
      sbg:fileTypes: TXT
      type: ['null', File]
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: 1, ramMin: 1000}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing:
      - {entry: "import argparse\n\nargs = argparse.ArgumentParser(description='This
          tool takes reference dictionary file as an input'\n                                             '
          and creates files which contain chromosome names grouped'\n                                             '
          based on their sizes.')\n\nargs.add_argument('--ref_dict', help='Reference
          dictionary', required=True)\nparsed = args.parse_args()\nref_dict = parsed.ref_dict\n\nwith
          open(ref_dict, 'r') as ref_dict_file:\n    sequence_tuple_list = []\n    longest_sequence
          = 0\n    for line in ref_dict_file:\n        if line.startswith(\"@SQ\"):\n
          \           line_split = line.split(\"\\t\")\n            # (Sequence_Name,
          Sequence_Length)\n            sequence_tuple_list.append((line_split[1].split(\"SN:\")[1],
          int(line_split[2].split(\"LN:\")[1])))\n    longest_sequence = sorted(sequence_tuple_list,
          key=lambda x: x[1], reverse=True)[0][1]\n# We are adding this to the intervals
          because hg38 has contigs named with embedded colons and a bug in GATK strips
          off\n# the last element after a :, so we add this as a sacrificial element.\nhg38_protection_tag
          = \":1+\"\n# initialize the tsv string with the first sequence\ntsv_string
          = sequence_tuple_list[0][0] + hg38_protection_tag\ntemp_size = sequence_tuple_list[0][1]\nfor
          sequence_tuple in sequence_tuple_list[1:]:\n    if temp_size + sequence_tuple[1]
          <= longest_sequence:\n        temp_size += sequence_tuple[1]\n        tsv_string
          += \"\\t\" + sequence_tuple[0] + hg38_protection_tag\n    else:\n        tsv_string
          += \"\\n\" + sequence_tuple[0] + hg38_protection_tag\n        temp_size
          = sequence_tuple[1]\n# add the unmapped sequences as a separate line to
          ensure that they are recalibrated as well\nwith open(\"./sequence_grouping.txt\",
          \"w\") as tsv_file:\n    tsv_file.write(tsv_string)\n    tsv_file.close()\n\ntsv_string
          += '\\n' + \"unmapped\"\n\nwith open(\"./sequence_grouping_with_unmapped.txt\",
          \"w\") as tsv_file_with_unmapped:\n    tsv_file_with_unmapped.write(tsv_string)\n
          \   tsv_file_with_unmapped.close()", entryname: CreateSequenceGroupingTSV.py,
        writable: false}
    - class: InlineJavascriptRequirement
      expressionLib: ["\nvar setMetadata = function(file, metadata) {\n    if (!('metadata'
          in file))\n        file['metadata'] = metadata;\n    else {\n        for
          (var key in metadata) {\n            file['metadata'][key] = metadata[key];\n
          \       }\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1,
          o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2
          = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example
          = o2[i]['metadata'];\n        for (var key in example) {\n            if
          (i == 0)\n                commonMetadata[key] = example[key];\n            else
          {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete
          commonMetadata[key]\n                }\n            }\n        }\n    }\n
          \   if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n
          \   } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i]
          = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BED Processing]
    sbg:content_hash: a9afa170a339934c60906ff616a6f2155426a9df80067bfc64f4140593aeffda6
    sbg:contributors: [nens, uros_sipetic]
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/6
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1555580154
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-createsequencegroupingtsv-4-1-0-0/4
    sbg:image_url: null
    sbg:latestRevision: 4
    sbg:license: BSD 3-clause
    sbg:links:
    - {id: 'https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels',
      label: GATK Germline GitHub}
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1558351560
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 4
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/6
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1555580154, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/1}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557734537, 'sbg:revision': 1, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/3}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557914517, 'sbg:revision': 2, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/4}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000609, 'sbg:revision': 3, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/5}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351560, 'sbg:revision': 4, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/createsequencegroupingtsv/6}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: 791.3856201171875
  sbg:y: 37.53794479370117
- id: gatk_gatherbamfiles_4_1_0_0
  in:
  - {default: true, id: create_index}
  - {id: output_prefix, source: output_prefix}
  - id: in_alignments
    source: [gatk_applybqsr_4_1_0_0/out_alignments]
  - {default: true, id: create_md5_file}
  label: GATK GatherBamFiles
  out:
  - {id: out_alignments}
  - {id: out_md5}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, shellQuote: false, valueFrom: /opt/gatk --java-options}
    - {position: 2, shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)
        {\n        return '\\\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\\\"';\n
        \   }\n    return '\\\"-Xmx2048M\\\"';\n}"}
    - {position: 4, shellQuote: false, valueFrom: "${\n    var tmp = [].concat(inputs.in_alignments);\n
        \       \n    if (inputs.output_prefix) {\n        return '-O ' +  inputs.output_prefix
        + \".bam\";\n        \n    }else if (tmp[0].metadata && tmp[0].metadata.sample_id)
        {\n        \n        return '-O ' +  tmp[0].metadata.sample_id + \".bam\";\n
        \   } else {\n         \n        return '-O ' +  tmp[0].path.split('/').pop().split(\".\")[0]
        + \".bam\";\n    }\n    \n    \n}"}
    - {position: 3, shellQuote: false, valueFrom: GatherBamFiles}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "**GATK GatherBamFiles** concatenates one or more BAM files resulted form
      scattered paralel anaysis. \n\n\n### Common Use Cases \n\n* **GATK GatherBamFiles**
      \ tool performs a rapid \"gather\" or concatenation on BAM files into single
      BAM file. This is often needed in operations that have been run in parallel
      across genomics regions by scattering their execution across computing nodes
      and cores thus resulting in smaller BAM files.\n* Usage example:\n```\n\njava
      -jar picard.jar GatherBamFiles\n      --INPUT=input1.bam\n      --INPUT=input2.bam\n```\n\n###
      Common Issues and Important Notes\n* **GATK GatherBamFiles** assumes that the
      list of BAM files provided as input are in the order that they should be concatenated
      and simply links the bodies of the BAM files while retaining the header from
      the first file. \n*  Operates by copying the gzip blocks directly for speed
      but also supports the generation of an MD5 in the output file and the indexing
      of the output BAM file.\n* This tool only support BAM files. It does not support
      SAM files.\n\n###Changes Intorduced by Seven Bridges\n* Generated output BAM
      file will be prefixed using the **Output prefix** parameter. In case the **Output
      prefix** is not provided, the output prefix will be the same as the **Sample
      ID** metadata from the **Input alignments**, if the **Sample ID** metadata exists.
      Otherwise, the output prefix will be inferred from the **Input alignments**
      filename. This way, having identical names of the output files between runs
      is avoided."
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-gatherbamfiles-4-1-0-0/9
    inputs:
    - doc: Memory overhead which will be allocated for one job.
      id: memory_overhead_per_job
      label: Memory Overhead Per Job
      sbg:category: Execution
      type: ['null', int]
    - doc: When writing files that need to be sorted, this will specify the number
        of records stored in ram before spilling to disk. Increasing this number reduces
        the number of file handles needed to sort the file, and increases the amount
        of ram needed.
      id: max_records_in_ram
      inputBinding: {position: 20, prefix: --MAX_RECORDS_IN_RAM, shellQuote: false}
      label: Max records in ram
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
      type: ['null', int]
    - doc: Memory which will be allocated for execution.
      id: memory_per_job
      label: Memory Per Job
      sbg:category: Execution
      type: ['null', int]
    - doc: Whether to create a bam index when writing a coordinate-sorted bam file.
      id: create_index
      inputBinding: {position: 4, prefix: --CREATE_INDEX, shellQuote: false}
      label: Create index
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Reference sequence file.
      id: in_reference
      inputBinding: {position: 7, prefix: --REFERENCE_SEQUENCE, shellQuote: false}
      label: Reference sequence
      sbg:altPrefix: -R
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', File]
    - doc: Name of the output bam file to write to.
      id: output_prefix
      label: Output prefix
      sbg:category: Optional Arguments
      type: ['null', string]
    - doc: Two or more bam files or text files containing lists of bam files (one
        per line). This argument must be specified at least once.
      id: in_alignments
      inputBinding: {position: 3, shellQuote: false, valueFrom: "${\n   if (self)\n
          \  {\n       var cmd = [];\n       for (var i = 0; i < self.length; i++)\n
          \      {\n           cmd.push('--INPUT', self[i].path);\n       }\n       return
          cmd.join(' ');\n   }\n\n}"}
      label: Input alignments
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM
      type: {items: File, type: array}
    - doc: Compression level for all compressed files created (e.g. Bam and vcf).
      id: compression_level
      inputBinding: {position: 4, prefix: --COMPRESSION_LEVEL, shellQuote: false}
      label: Compression level
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
      type: ['null', int]
    - doc: Validation stringency for all sam files read by this program. Setting stringency
        to silent can improve performance when processing a bam file in which variable-length
        data (read, qualities, tags) do not otherwise need to be decoded.
      id: validation_stringency
      inputBinding: {position: 4, prefix: --VALIDATION_STRINGENCY, shellQuote: false}
      label: Validation stringency
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
      type:
      - 'null'
      - name: validation_stringency
        symbols: [STRICT, LENIENT, SILENT]
        type: enum
    - doc: Whether to create an MD5 digest for any BAM or FASTQ files created.
      id: create_md5_file
      inputBinding: {position: 5, prefix: --CREATE_MD5_FILE, shellQuote: false}
      label: Create MD5 file
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'FALSE'
      type: ['null', boolean]
    - doc: This input allows a user to set the desired CPU requirement when running
        a tool or adding it to a workflow.
      id: cpu_per_job
      label: CPU per job
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    label: GATK GatherBamFiles
    outputs:
    - doc: Output BAM file obtained by merging input BAM files.
      id: out_alignments
      label: Output BAM file
      outputBinding: {glob: '*.bam', outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: BAM
      secondaryFiles: ["${\n    if (inputs.create_index)\n    {\n        return [self.basename
          + \".bai\", self.nameroot + \".bai\"];\n    }\n    else {\n        return
          ''; \n    }\n}"]
      type: ['null', File]
    - doc: MD5 ouput BAM file.
      id: out_md5
      label: MD5 file
      outputBinding: {glob: '*.md5', outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: MD5
      type: ['null', File]
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: "${\n    return inputs.cpu_per_job ?
        inputs.cpu_per_job : 1;\n}", ramMin: "${\n    var memory = 4096;\n    if (inputs.memory_per_job)
        \n    {\n        memory = inputs.memory_per_job;\n    }\n    if (inputs.memory_overhead_per_job)\n
        \   {\n        memory += inputs.memory_overhead_per_job;\n    }\n    return
        memory;\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n
          \   for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n
          \   }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n
          \   var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2
          = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example
          = o2[i]['metadata'];\n        for (var key in example) {\n            if
          (i == 0)\n                commonMetadata[key] = example[key];\n            else
          {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete
          commonMetadata[key]\n                }\n            }\n        }\n    }\n
          \   if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n
          \   } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i]
          = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar
          toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy
          = function(files, key) {\n    var groupedFiles = [];\n    var tempDict =
          {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n
          \       if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};", "\nvar setMetadata = function(file,
          metadata) {\n    if (!('metadata' in file))\n        file['metadata'] =
          metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BAM Processing]
    sbg:content_hash: adc3fdd806bf7e70cfd29e650f70e8bdc6477baa1d0dc7ef7792f2f8806bcd064
    sbg:contributors: [nens]
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/23
    sbg:createdBy: nens
    sbg:createdOn: 1554894822
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-gatherbamfiles-4-1-0-0/9
    sbg:image_url: null
    sbg:latestRevision: 9
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - {id: 'https://software.broadinstitute.org/gatk/', label: Homepage}
    - {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_GatherBamFiles.php',
      label: Documentation}
    - {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
    - {id: 'https://github.com/broadinstitute/gatk/', label: Source}
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1558531990
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 9
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/23
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1554894822, 'sbg:revision': 0, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/11}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557734548, 'sbg:revision': 1, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/14}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557914509, 'sbg:revision': 2, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/16}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000604, 'sbg:revision': 3, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/17}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351555, 'sbg:revision': 4, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/18}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558451620, 'sbg:revision': 5, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/19}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558525775, 'sbg:revision': 6, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/20}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558526183, 'sbg:revision': 7, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/21}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558528334, 'sbg:revision': 8, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/22}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558531990, 'sbg:revision': 9, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbamfiles-4-1-0-0/23}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: 1867.5662841796875
  sbg:y: 208.6806640625
- id: gatk_gatherbqsrreports_4_1_0_0
  in:
  - id: in_bqsr_reports
    source: [gatk_baserecalibrator_4_1_0_0/output]
  label: GATK GatherBQSRReports
  out:
  - {id: out_gathered_bqsr_reports}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, shellQuote: false, valueFrom: /opt/gatk}
    - {position: 1, shellQuote: false, valueFrom: --java-options}
    - {position: 2, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)
        {\n        return '\\\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\\\"';\n
        \   }\n    return '\\\"-Xmx2048M\\\"';\n}"}
    - {position: 3, shellQuote: false, valueFrom: GatherBQSRReports}
    - {position: 6, prefix: '', shellQuote: false, valueFrom: "${\n    var tmp = [].concat(inputs.in_bqsr_reports);\n
        \       \n    if (inputs.output_prefix) {\n        return '-O ' +  inputs.output_prefix
        + \".recal_data.csv\";\n        \n    }else if (tmp[0].metadata && tmp[0].metadata.sample_id)
        {\n        \n        return '-O ' +  tmp[0].metadata.sample_id + \".recal_data.csv\";\n
        \   } else {\n         \n        return '-O ' +  tmp[0].path.split('/').pop().split(\".\")[0]
        + \".recal_data.csv\";\n    }\n    \n    \n}"}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "**GATK GatherBQSRReports** gathers scattered BQSR recalibration reports
      from parallelized base recalibration runs into a single file.\nThe combination
      is done simply by adding up all observations and errors.\n\n\n### Common Use
      Cases \n\n* This tool is intended to be used to combine recalibration tables
      from runs of BaseRecalibrator parallelized per-interval.\n\n* Usage example
      (input BAM file is single-end):\n\n```\ngatk GAtherBQSRReports \n     --input
      example1.csv\n     --input example2.csv\n```\n\n\n### Common Issues and Important
      Notes\n\n* This method DOES NOT recalculate the empirical qualities and quantized
      qualities. You have to recalculate them after combining. The reason for not
      calculating it is because this function is intended for combining a series of
      recalibration reports, and it only makes sense to calculate the empirical qualities
      and quantized qualities after all the recalibration reports have been combined.
      This is done to make the tool faster.\nThe reported empirical quality is recalculated.\n\n###Changes
      Introduced by Seven Bridges\n\n* Output file will be prefixed using the **Output
      prefix** parameter. In case **Output prefix** is not provided, output prefix
      will be the same as the **Sample ID** metadata from the **Input bqsr reports**,
      if the **Sample ID** metadata exists. Otherwise, output prefix will be inferred
      from the **Input bqsr reports** filename. This way, having identical names of
      the output files between runs is avoided. Moreover,  **.recal_data.csv** will
      be added before the extension of the output file name. \n\n\n###References\n\n[1]
      [GATK GatherBQSRReports](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_GatherBQSRReports.php)"
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-gatherbqsrreports-4-1-0-0/7
    inputs:
    - doc: List of scattered BQSR report files. This argument must be specified at
        least once.
      id: in_bqsr_reports
      inputBinding: {itemSeparator: 'null', position: 4, prefix: '', shellQuote: false,
        valueFrom: "${\n   if (self)\n   {\n       var cmd = [];\n       for (var
          i = 0; i < self.length; i++)\n       {\n           cmd.push('--input', self[i].path);\n
          \      }\n       return cmd.join(' ');\n   }\n\n}"}
      label: Input bqsr reports
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: CSV
      type: {items: File, type: array}
    - doc: Output prefix for the  gathered bqsr report.
      id: output_prefix
      label: Gathered report output prefix
      sbg:category: Optional argument
      sbg:toolDefaultValue: Inferred from sample metadata.
      type: ['null', string]
    - doc: Memory which will be allocated for execution.
      id: memory_per_job
      label: Memory Per Job
      sbg:category: Execution
      type: ['null', int]
    - doc: Memory overhead which will be allocated for one job.
      id: memory_overhead_per_job
      label: Memory Overhead Per Job
      sbg:category: Execution
      type: ['null', int]
    - doc: Number of CPUs which will be allocated for the job.
      id: cpu_per_job
      label: CPU per job
      sbg:category: Execution
      type: ['null', int]
    label: GATK GatherBQSRReports
    outputs:
    - doc: Gathered report - merged inputs.
      id: out_gathered_bqsr_reports
      label: Gathered report
      outputBinding: {glob: '*.csv', outputEval: '$(inheritMetadata(self, inputs.in_bqsr_reports))'}
      sbg:fileTypes: CSV
      type: ['null', File]
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: '$(inputs.cpu_per_job ? inputs.cpu_per_job
        : 1)', ramMin: "${\n    if (inputs.memory_per_job) {\n        if (inputs.memory_overhead_per_job)
        {\n            return inputs.memory_per_job + inputs.memory_overhead_per_job;\n
        \       } \n        else {\n            return inputs.memory_per_job;\n        }\n
        \   } \n    \n    else if (!inputs.memory_per_job && inputs.memory_overhead_per_job)
        {\n        return 2048 + inputs.memory_overhead_per_job;\n    } \n    \n    else
        {\n        return 2048;\n    }\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n
          \   for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n
          \   }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n
          \   var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2
          = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example
          = o2[i]['metadata'];\n        for (var key in example) {\n            if
          (i == 0)\n                commonMetadata[key] = example[key];\n            else
          {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete
          commonMetadata[key]\n                }\n            }\n        }\n    }\n
          \   if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n
          \   } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i]
          = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar
          toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy
          = function(files, key) {\n    var groupedFiles = [];\n    var tempDict =
          {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n
          \       if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};", "\nvar setMetadata = function(file,
          metadata) {\n    if (!('metadata' in file))\n        file['metadata'] =
          metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BAM Processing]
    sbg:content_hash: a4a6b8fb1bbf940a7e1d86e48fedb7300099bd07c50720bd9d07c55192a893565
    sbg:contributors: [nens, uros_sipetic]
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/25
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1554810073
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-gatherbqsrreports-4-1-0-0/7
    sbg:image_url: null
    sbg:latestRevision: 7
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - {id: 'https://software.broadinstitute.org/gatk/', label: Homepage}
    - {id: 'https://github.com/broadinstitute/gatk/', label: Source code}
    - {id: 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip',
      label: Download}
    - {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
    - {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_GatherBQSRReports.php',
      label: Documentation}
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1558451160
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 7
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/25
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1554810073, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/8}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1554894740, 'sbg:revision': 1, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/11}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557487015, 'sbg:revision': 2, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/13}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557734524, 'sbg:revision': 3, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/17}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557744219, 'sbg:revision': 4, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/22}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000599, 'sbg:revision': 5, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/23}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351550, 'sbg:revision': 6, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/24}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558451160, 'sbg:revision': 7, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-gatherbqsrreports-4-1-0-0/25}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: 1459.721923828125
  sbg:y: 330.2196960449219
- id: gatk_markduplicates_4_1_0_0
  in:
  - {default: queryname, id: assume_sort_order}
  - id: in_alignments
    source: [gatk_mergebamalignment_4_1_0_0/out_alignments]
  - {default: 2500, id: optical_duplicate_pixel_distance}
  - {default: SILENT, id: validation_stringency}
  label: GATK MarkDuplicates
  out:
  - {id: out_alignments}
  - {id: output_metrics}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, prefix: '', shellQuote: false, valueFrom: /opt/gatk}
    - {position: 1, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)\n
        \   {\n        return \"--java-options\";\n    }\n    else {\n        return
        ''; \n    }\n}\n    "}
    - {position: 2, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)
        {\n        return '\\\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\\\"';\n
        \   }\n    else {\n        return ''; \n    }\n}"}
    - {position: 3, shellQuote: false, valueFrom: MarkDuplicates}
    - {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n    var in_alignments
        = [].concat(inputs.in_alignments);\n    var output_ext = inputs.output_file_format
        ? \".\" + inputs.output_file_format : in_alignments[0].nameext;\n    var output_prefix
        = '';\n    if (inputs.output_prefix)\n    {\n        output_prefix = inputs.output_prefix;\n
        \   }\n    else\n    {\n        if (in_alignments[0].metadata && in_alignments[0].metadata.sample_id)\n
        \       {\n            output_prefix = in_alignments[0].metadata.sample_id;\n
        \       }\n        else\n        {\n            output_prefix = in_alignments[0].nameroot.split('.')[0];\n
        \       }\n    }\n    return \"--OUTPUT \" + output_prefix + \".dedupped\"
        + output_ext;\n}"}
    - {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n    var in_alignments
        = [].concat(inputs.in_alignments);\n    var output_prefix = '';  \n\n    if
        (inputs.output_prefix)\n    {\n        output_prefix = inputs.output_prefix;\n
        \   }\n    else\n    {\n        if (in_alignments[0].metadata && in_alignments[0].metadata.sample_id)\n
        \       {\n            output_prefix = in_alignments[0].metadata.sample_id;\n
        \       }\n        else\n        {\n            output_prefix = in_alignments[0].nameroot.split('.')[0];\n
        \       }\n    }\n    return \"--METRICS_FILE \" + output_prefix + \".dedupped.metrics\";\n}"}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "The **GATK  MarkDuplicates** tool identifies duplicate reads in a BAM or
      SAM file.\n\nThis tool locates and tags duplicate reads in a BAM or SAM file,
      where duplicate reads are defined as originating from a single fragment of DNA.
      Duplicates can arise during sample preparation e.g. library construction using
      PCR. Duplicate reads can also result from a single amplification cluster, incorrectly
      detected as multiple clusters by the optical sensor of the sequencing instrument.
      These duplication artifacts are referred to as optical duplicates [1].\n\nThe
      MarkDuplicates tool works by comparing sequences in the 5 prime positions of
      both reads and read-pairs in the SAM/BAM file. The **Barcode tag** (`--BARCODE_TAG`)
      option is available to facilitate duplicate marking using molecular barcodes.
      After duplicate reads are collected, the tool differentiates the primary and
      duplicate reads using an algorithm that ranks reads by the sums of their base-quality
      scores (default method).\n\n\n###Common Use Cases\n\n* The **GATK MarkDuplicates**
      tool requires the BAM or SAM file on its **Input BAM/SAM file** (`--INPUT`)
      input. The tool generates a new SAM or BAM file on its **Output BAM/SAM** output,
      in which duplicates have been identified in the SAM flags field for each read.
      Duplicates are marked with the hexadecimal value of 0x0400, which corresponds
      to a decimal value of 1024. If you are not familiar with this type of annotation,
      please see the following [blog post](https://software.broadinstitute.org/gatk/blog?id=7019)
      for additional information. **MarkDuplicates** also produces a metrics file
      on its **Output metrics file** output, indicating the numbers of duplicates
      for both single and paired end reads.\n\n* The program can take either coordinate-sorted
      or query-sorted inputs, however the behavior is slightly different. When the
      input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary
      alignments are not marked as duplicates. However, when the input is query-sorted
      (actually query-grouped), then unmapped mates and secondary/supplementary reads
      are not excluded from the duplication test and can be marked as duplicate reads.\n\n*
      If desired, duplicates can be removed using the **Remove duplicates** (`--REMOVE_DUPLICATES`)
      and **Remove sequencing duplicates** ( `--REMOVE_SEQUENCING_DUPLICATES`) options.\n\n*
      Although the bitwise flag annotation indicates whether a read was marked as
      a duplicate, it does not identify the type of duplicate. To do this, a new tag
      called the duplicate type (DT) tag was recently added as an optional output
      of a SAM/BAM file. Invoking the **Tagging policy** ( `--TAGGING_POLICY`) option,
      you can instruct the program to mark all the duplicates (All), only the optical
      duplicates (OpticalOnly), or no duplicates (DontTag). The records within the
      output SAM/BAM file will have values for the 'DT' tag (depending on the invoked
      **TAGGING_POLICY** option), as either library/PCR-generated duplicates (LB),
      or sequencing-platform artifact duplicates (SQ). \n\n* This tool uses the **Read
      name regex** (`--READ_NAME_REGEX`) and the **Optical duplicate pixel distance**
      (`--OPTICAL_DUPLICATE_PIXEL_DISTANCE`) options as the primary methods to identify
      and differentiate duplicate types. Set **READ_NAME_REGEX** to null to skip optical
      duplicate detection, e.g. for RNA-seq or other data where duplicate sets are
      extremely large and estimating library complexity is not an aim. Note that without
      optical duplicate counts, library size estimation will be inaccurate.\n\n* Usage
      example:\n\n```\ngatk MarkDuplicates \\\n      --INPUT input.bam \\\n      --OUTPUT
      marked_duplicates.bam \\\n      --METRICS_FILE marked_dup_metrics.txt\n```\n\n###Changes
      Introduced by Seven Bridges\n\n* All output files will be prefixed using the
      **Output prefix** parameter. In case **Output prefix** is not provided, output
      prefix will be the same as the Sample ID metadata from the **Input SAM/BAM file**,
      if the Sample ID metadata exists. Otherwise, output prefix will be inferred
      from the **Input SAM/BAM** filename. This way, having identical names of the
      output files between runs is avoided. Moreover,  **dedupped** will be added
      before the extension of the output file name. \n\n* The user has a possibility
      to specify the output file format using the **Output file format** option. Otherwise,
      the output file format will be the same as the format of the input file.\n\n###Common
      Issues and Important Notes\n\n* None\n\n###Performance Benchmarking\n\nBelow
      is a table describing runtimes and task costs of **GATK MarkDuplicates** for
      a couple of different samples, executed on the AWS cloud instances:\n\n| Experiment
      type |  Input size | Duration |  Cost | Instance (AWS) | \n|:--------------:|:------------:|:--------:|:-------:|:---------:|\n|
      \    RNA-Seq     |  1.8 GB |   3min   | ~0.02$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq
      \    |  5.3 GB |   9min   | ~0.06$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq     |
      8.8 GB |  16min  | ~0.11$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq     | 17 GB
      |  30min  | ~0.20$ | c4.2xlarge (8 CPUs) |\n\n*Cost can be significantly reduced
      by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances)
      for more details.*\n\n###References\n\n[1] [GATK MarkDuplicates](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_markduplicates_MarkDuplicates.php)"
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-markduplicates-4-1-0-0/12
    inputs:
    - doc: Add PG tag to each read in a SAM or BAM file.
      id: add_pg_tag_to_reads
      inputBinding: {position: 4, prefix: --ADD_PG_TAG_TO_READS, shellQuote: false}
      label: Add PG tag to reads
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: add_pg_tag_to_reads
        symbols: ['true', 'false']
        type: enum
    - doc: If not null, assume that the input file has this order even if the header
        says otherwise. Cannot be used in conjuction with argument(s) ASSUME_SORTED
        (AS).
      id: assume_sort_order
      inputBinding: {position: 4, prefix: --ASSUME_SORT_ORDER, shellQuote: false}
      label: Assume sort order
      sbg:altPrefix: -ASO
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - name: assume_sort_order
        symbols: [unsorted, queryname, coordinate, duplicate, unknown]
        type: enum
    - doc: 'If true, assume that the input file is coordinate sorted even if the header
        says otherwise. Deprecated, used ASSUME_SORT_ORDER=coordinate instead. Exclusion:
        This argument cannot be used at the same time as ASSUME_SORT_ORDER (ASO).'
      id: assume_sorted
      inputBinding: {position: 4, prefix: --ASSUME_SORTED, shellQuote: false}
      label: Assume sorted
      sbg:altPrefix: -AS
      sbg:category: Optional arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Barcode SAM tag (ex. BC for 10x genomics).
      id: barcode_tag
      inputBinding: {position: 4, prefix: --BARCODE_TAG, shellQuote: false}
      label: Barcode tag
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: Clear DT tag from input SAM records. Should be set to false if input SAM
        doesn't have this tag.
      id: clear_dt
      inputBinding: {position: 4, prefix: --CLEAR_DT, shellQuote: false}
      label: Clear DT
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: clear_dt
        symbols: ['true', 'false']
        type: enum
    - doc: Comment(s) to include in the output file's header.
      id: comment
      inputBinding: {position: 4, shellQuote: false, valueFrom: "${\n    if (self)\n
          \   {\n        var cmd = [];\n        for (var i = 0; i < self.length; i++)
          \n        {\n            cmd.push('--COMMENT', self[i]);\n            \n
          \       }\n        return cmd.join(' ');\n    }\n}"}
      label: Comment
      sbg:altPrefix: -CO
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Compression level for all compressed files created (e.g. BAM and VCF).
      id: compression_level
      inputBinding: {position: 4, prefix: --COMPRESSION_LEVEL, shellQuote: false}
      label: Compression level
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
      type: ['null', int]
    - doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
      id: create_index
      inputBinding: {position: 4, prefix: --CREATE_INDEX, shellQuote: false}
      label: Create index
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Treat UMIs as being duplex stranded. This option requires that the UMI
        consist of two equal length strings that are separated by a hyphen (e.g. 'ATC-GTC').
        Reads are considered duplicates if, in addition to standard definition, have
        identical normalized UMIs. A UMI from the 'bottom' strand is normalized by
        swapping its content around the hyphen (eg. ATC-GTC becomes GTC-ATC). A UMI
        from the 'top' strand is already normalized as it is. Both reads from a read
        pair considered top strand if the read 1 unclipped 5' coordinate is less than
        the read 2 unclipped 5' coordinate. All chimeric reads and read fragments
        are treated as having come from the top strand. With this option it is required
        that the BARCODE_TAG hold non-normalized UMIs.
      id: duplex_umi
      inputBinding: {position: 4, prefix: --DUPLEX_UMI, shellQuote: false}
      label: Duplex UMI
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: The scoring strategy for choosing the non-duplicate among candidates.
      id: duplicate_scoring_strategy
      inputBinding: {position: 4, prefix: --DUPLICATE_SCORING_STRATEGY, shellQuote: false}
      label: Duplicate scoring strategy
      sbg:altPrefix: -DS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: SUM_OF_BASE_QUALITIES
      type:
      - 'null'
      - name: duplicate_scoring_strategy
        symbols: [SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM]
        type: enum
    - doc: Input SAM or BAM files to analyze. Must be coordinate sorted.
      id: in_alignments
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   var in_files = [].concat(inputs.in_alignments);\n    if (in_files)\n
          \   {\n        var cmd = [];\n        for (var i = 0; i < in_files.length;
          i++) \n        {\n            cmd.push('--INPUT', in_files[i].path);\n        }\n
          \       return cmd.join(' ');\n    }\n}"}
      label: Input BAM/SAM file
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM
      type: {items: File, type: array}
    - doc: Maximum number of file handles to keep open when spilling read ends to
        disk. Set this number a little lower than the per-process maximum number of
        file that may be open. This number can be found by executing the 'ulimit -n'
        command on a unix system.
      id: max_file_handles_for_read_ends_map
      inputBinding: {position: 4, prefix: --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP, shellQuote: false}
      label: Max file handles for read ends map
      sbg:altPrefix: -MAX_FILE_HANDLES
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '8000'
      type: ['null', int]
    - doc: This number is the maximum size of a set of duplicate reads for which we
        will attempt to determine which are optical duplicates. Please be aware that
        if you raise this value too high and do encounter a very large set of duplicate
        reads, it will severely affect the runtime of this tool. To completely disable
        this check, set the value to -1.
      id: max_optical_duplicate_set_size
      inputBinding: {position: 4, prefix: --MAX_OPTICAL_DUPLICATE_SET_SIZE, shellQuote: false}
      label: Max optical duplicate set size
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '300000'
      type: ['null', int]
    - doc: When writing files that need to be sorted, this will specify the number
        of records stored in RAM before spilling to disk. Increasing this number reduces
        the number of file handles needed to sort the file, and increases the amount
        of RAM needed.
      id: max_records_in_ram
      inputBinding: {position: 4, prefix: --MAX_RECORDS_IN_RAM, shellQuote: false}
      label: Max records in RAM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
      type: ['null', int]
    - doc: This input allows a user to set the desired overhead memory when running
        a tool or adding it to a workflow. This amount will be added to the Memory
        per job in the Memory requirements section but it will not be added to the
        -Xmx parameter leaving some memory not occupied which can be used as stack
        memory (-Xmx parameter defines heap memory). This input should be defined
        in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      id: memory_overhead_per_job
      label: Memory overhead per job
      sbg:category: Platform Options
      type: ['null', int]
    - doc: This input allows a user to set the desired memory requirement when running
        a tool or adding it to a workflow. This value should be propagated to the
        -Xmx parameter too.This input should be defined in MB (for both the platform
        part and the -Xmx part if Java tool is wrapped).
      id: memory_per_job
      label: Memory per job
      sbg:category: Platform Options
      type: ['null', int]
    - doc: SAM tag to uniquely identify the molecule from which a read was derived.
        Use of this option requires that the BARCODE_TAG option be set to a non null
        value.
      id: molecular_identifier_tag
      inputBinding: {position: 4, prefix: --MOLECULAR_IDENTIFIER_TAG, shellQuote: false}
      label: Molecular identifier tag
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: The maximum offset between two duplicate clusters in order to consider
        them optical duplicates. The default is appropriate for unpatterned versions
        of the illumina platform. For the patterned flowcell models, 2500 is moreappropriate.
        For other platforms and models, users should experiment to find what works
        best.
      id: optical_duplicate_pixel_distance
      inputBinding: {position: 4, prefix: --OPTICAL_DUPLICATE_PIXEL_DISTANCE, shellQuote: false}
      label: Optical duplicate pixel distance
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '100'
      type: ['null', int]
    - doc: Value of CL tag of PG record to be created. If not supplied the command
        line will be detected automatically.
      id: program_group_command_line
      inputBinding: {position: 4, prefix: --PROGRAM_GROUP_COMMAND_LINE, shellQuote: false}
      label: Program group command line
      sbg:altPrefix: -PG_COMMAND
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: Value of PN tag of PG record to be created.
      id: program_group_name
      inputBinding: {position: 4, prefix: --PROGRAM_GROUP_NAME, shellQuote: false}
      label: Program group name
      sbg:altPrefix: -PG_NAME
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: MarkDuplicates
      type: ['null', string]
    - doc: Value of VN tag of PG record to be created. If not specified, the version
        will be detected automatically.
      id: program_group_version
      inputBinding: {position: 4, prefix: --PROGRAM_GROUP_VERSION, shellQuote: false}
      label: Program group version
      sbg:altPrefix: -PG_VERSION
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: The program record ID for the @PG record(s) created by this program. Set
        to null to disable PG record creation.  This string may have a suffix appended
        to avoid collision with other program record IDs.
      id: program_record_id
      inputBinding: {position: 4, prefix: --PROGRAM_RECORD_ID, shellQuote: false}
      label: Program record id
      sbg:altPrefix: -PG
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: MarkDuplicates
      type: ['null', string]
    - doc: 'MarkDuplicates can use the tile and cluster positions to estimate the
        rate of optical duplication in addition to the dominant source of duplication,
        PCR, to provide a more accurate estimation of library size. By default (with
        no READ_NAME_REGEX specified), MarkDuplicates will attempt to extract coordinates
        using a split on '':'' (see note below). Set READ_NAME_REGEX to ''null'' to
        disable optical duplicate detection. Note that without optical duplicate counts,
        library size estimation will be less accurate. If the read name does not follow
        a standard illumina colon-separation convention, but does contain tile and
        x,y coordinates, a regular expression can be specified to extract three variables:
        tile/region, x coordinate and y coordinate from a read name. The regular expression
        must contain three capture groups for the three variables, in order. It must
        match the entire read name. e.g. if field names were separated by semi-colon
        ('';'') this example regex could be specified (?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$
        Note that if no READ_NAME_REGEX is specified, the read name is split on '':''.
        For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile,
        x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements
        are assumed to be tile, x and y values.'
      id: read_name_regex
      inputBinding: {position: 4, prefix: --READ_NAME_REGEX, shellQuote: false}
      label: Read name regex
      sbg:category: Optional Arguments
      type: ['null', string]
    - doc: Read one barcode SAM tag (ex. BX for 10x Genomics).
      id: read_one_barcode_tag
      inputBinding: {position: 4, prefix: --READ_ONE_BARCODE_TAG, shellQuote: false}
      label: Read one barcode tag
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: Read two barcode SAM tag (ex. BX for 10x Genomics).
      id: read_two_barcode_tag
      inputBinding: {position: 4, prefix: --READ_TWO_BARCODE_TAG, shellQuote: false}
      label: Read two barcode tag
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: If true do not write duplicates to the output file instead of writing them
        with appropriate flags set.
      id: remove_duplicates
      inputBinding: {position: 4, prefix: --REMOVE_DUPLICATES, shellQuote: false}
      label: Remove duplicates
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: If true remove 'optical' duplicates and other duplicates that appear to
        have arisen from the sequencing process instead of the library preparation
        process, even if REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true,
        all duplicates are removed and this option is ignored.
      id: remove_sequencing_duplicates
      inputBinding: {position: 4, prefix: --REMOVE_SEQUENCING_DUPLICATES, shellQuote: false}
      label: Remove sequencing duplicates
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: This number, plus the maximum RAM available to the JVM, determine the memory
        footprint used by some of the sorting collections. If you are running out
        of memory, try reducing this number.
      id: sorting_collection_size_ratio
      inputBinding: {position: 4, prefix: --SORTING_COLLECTION_SIZE_RATIO, shellQuote: false}
      label: Sorting collection size ratio
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0.25'
      type: ['null', float]
    - doc: If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG
        (DS), indicates the size of the duplicate set. The smallest possible DS value
        is 2 which occurs when two reads map to the same portion of the reference
        only one of which is marked as duplicate. The second tag, DUPLICATE_SET_INDEX_TAG
        (DI), represents a unique identifier for the duplicate set to which the record
        belongs. This identifier is the index-in-file of the representative read that
        was selected out of the duplicate set.
      id: tag_duplicate_set_members
      inputBinding: {position: 4, prefix: --TAG_DUPLICATE_SET_MEMBERS, shellQuote: false}
      label: Tag duplicate set members
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Determines how duplicate types are recorded in the DT optional attribute.
      id: tagging_policy
      inputBinding: {position: 4, prefix: --TAGGING_POLICY, shellQuote: false}
      label: Tagging policy
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: DontTag
      type:
      - 'null'
      - name: tagging_policy
        symbols: [DontTag, OpticalOnly, All]
        type: enum
    - doc: Validation stringency for all SAM files read by this program. Setting stringency
        to SILENT can improve performance when processing a BAM file in which variable-length
        data (read, qualities, tags) do not otherwise need to be decoded.
      id: validation_stringency
      inputBinding: {position: 4, prefix: --VALIDATION_STRINGENCY, shellQuote: false}
      label: Validation stringency
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
      type:
      - 'null'
      - name: validation_stringency
        symbols: [STRICT, LENIENT, SILENT]
        type: enum
    - doc: Output file name prefix.
      id: output_prefix
      label: Output prefix
      sbg:category: Optional Arguments
      type: ['null', string]
    - doc: Output file format
      id: output_file_format
      label: Output file format
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: BAM
      type:
      - 'null'
      - name: output_file_format
        symbols: [bam, sam]
        type: enum
    - doc: This input allows a user to set the desired CPU requirement when running
        a tool or adding it to a workflow.
      id: cpu_per_job
      label: CPU per job
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    label: GATK MarkDuplicates
    outputs:
    - doc: Output BAM/SAM file which contains marked records.
      id: out_alignments
      label: Output BAM/SAM file
      outputBinding: {glob: '*am', outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: BAM, SAM
      secondaryFiles: ["${ \r\n   if (inputs.create_index)   {\r\n       return [self.basename
          + \".bai\", self.nameroot + \".bai\"]\r\n   }  else {\r\n       return [];
          \r\n  }\r\n}"]
      type: ['null', File]
    - doc: Output duplication metrics file.
      id: output_metrics
      label: Output metrics file
      outputBinding: {glob: '*metrics', outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: METRICS
      type: File
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: "${\n    return inputs.cpu_per_job ?
        inputs.cpu_per_job : 1;\n}", ramMin: "${\n    var memory = 4096;\n    if (inputs.memory_per_job)
        \n    {\n        memory = inputs.memory_per_job;\n    }\n    if (inputs.memory_overhead_per_job)\n
        \   {\n        memory += inputs.memory_overhead_per_job;\n    }\n    return
        memory;\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n
          \   else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar
          groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict
          = {};\n    for (var i = 0; i < files.length; i++) {\n        var value =
          files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};", "\nvar setMetadata = function(file,
          metadata) {\n    if (!('metadata' in file))\n        file['metadata'] =
          metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BAM Processing]
    sbg:content_hash: a112438cd40b078b2fbf816496a7cabec5688e19c781aac7f79a1de917e0eabfb
    sbg:contributors: [uros_sipetic, nemanja.vucic, veliborka_josipovic, nens]
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/26
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552668097
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-markduplicates-4-1-0-0/12
    sbg:image_url: null
    sbg:latestRevision: 12
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - {id: 'https://software.broadinstitute.org/gatk/', label: Homepage}
    - {id: 'https://github.com/broadinstitute/gatk/', label: Source Code}
    - {id: 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip',
      label: Download}
    - {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
    - {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_markduplicates_MarkDuplicates.php',
      label: Documentation}
    sbg:modifiedBy: uros_sipetic
    sbg:modifiedOn: 1562416183
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 12
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/26
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1552668097, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/9}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554492835, 'sbg:revision': 1,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/13}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554720881, 'sbg:revision': 2,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/14}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554999255, 'sbg:revision': 3,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/15}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1555945044, 'sbg:revision': 4,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/17}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557734534, 'sbg:revision': 5, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/18}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000580, 'sbg:revision': 6, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/19}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351536, 'sbg:revision': 7, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/21}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558447931, 'sbg:revision': 8, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/22}
    - {'sbg:modifiedBy': nemanja.vucic, 'sbg:modifiedOn': 1559750423, 'sbg:revision': 9,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/23}
    - {'sbg:modifiedBy': nemanja.vucic, 'sbg:modifiedOn': 1559751034, 'sbg:revision': 10,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/24}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1561632463, 'sbg:revision': 11, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/25}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1562416183, 'sbg:revision': 12,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-markduplicates-4-1-0-0/26}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: 255
  sbg:y: 102.40715789794922
- id: gatk_mergebamalignment_4_1_0_0
  in:
  - {default: 'true', id: add_mate_cigar}
  - id: in_alignments
    source: [bwa_mem_bundle_0_7_15/aligned_reads]
    valueFrom: $([self])
  - {default: true, id: aligner_proper_pair_flags}
  - default: [X0]
    id: attributes_to_retain
  - {default: 'false', id: clip_adapters}
  - default: [FR]
    id: expected_orientations
  - {default: -1, id: max_insertions_or_deletions}
  - {default: 2000000, id: max_records_in_ram, source: max_records_in_ram}
  - {default: 'true', id: paired_run}
  - {default: MostDistant, id: primary_alignment_strategy}
  - {default: '"bwa mem -K 100000000 -p -v 3 -t 16 -Y ref_fasta"', id: program_group_command_line}
  - {default: bwamem, id: program_group_name}
  - {default: 0.7.15, id: program_group_version}
  - {default: bwamem, id: program_record_id}
  - {id: in_reference, source: in_reference}
  - {default: unsorted, id: sort_order}
  - {default: true, id: unmap_contaminant_reads}
  - {id: unmapped_bam, source: in_alignments}
  - {default: COPY_TO_TAG, id: unmapped_read_strategy}
  - {default: SILENT, id: validation_stringency}
  label: GATK MergeBamAlignment
  out:
  - {id: out_alignments}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, prefix: '', shellQuote: false, valueFrom: /opt/gatk}
    - {position: 1, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)\n
        \   {\n        return \"--java-options\";\n    }\n    else {\n        return
        '';\n    }\n}"}
    - {position: 2, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)
        {\n        return '\\\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\\\"';\n
        \   }\n    else {\n        return ''; \n        \n    }\n}"}
    - {position: 3, shellQuote: false, valueFrom: MergeBamAlignment}
    - {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n    var in_alignments
        = [].concat(inputs.in_alignments);\n    var output_ext = inputs.output_file_format
        ? inputs.output_file_format : in_alignments[0].path.split('.').pop();\n    var
        output_prefix = '';\n    var file1_name = ''; \n    var file2_name = ''; \n
        \   if (inputs.output_prefix)\n    {\n        output_prefix = inputs.output_prefix;\n
        \   }\n    else \n    {\n        if (in_alignments.length > 1)\n        {\n
        \           in_alignments.sort(function(file1, file2) {\n                file1_name
        = file1.path.split('/').pop().toUpperCase();\n                file2_name =
        file2.path.split('/').pop().toUpperCase();\n                if (file1_name
        < file2_name) {\n                    return -1;\n                }\n                if
        (file1_name > file2_name) {\n                    return 1;\n                }\n
        \               // names must be equal\n                return 0;\n            });\n
        \       }\n        \n        var in_alignments_first =  in_alignments[0];\n
        \       if (in_alignments_first.metadata && in_alignments_first.metadata.sample_id)\n
        \       {\n            output_prefix = in_alignments_first.metadata.sample_id;\n
        \       }\n        else \n        {\n            output_prefix = in_alignments_first.path.split('/').pop().split('.')[0];\n
        \       }\n        \n        if (in_alignments.length > 1)\n        {\n            output_prefix
        = output_prefix + \".\" + in_alignments.length;\n        }\n    }\n    \n
        \   return \"--OUTPUT \" + output_prefix + \".merged.\" + output_ext;\n}"}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "The **GATK MergeBamAlignment** tool is used for merging BAM/SAM alignment
      info from a third-party aligner with the data in an unmapped BAM file, producing
      a third BAM file that has alignment data (from the aligner) and all the remaining
      data from the unmapped BAM.\n\nMany alignment tools still require FASTQ format
      input. The unmapped BAM may contain useful information that will be lost in
      the conversion to FASTQ (meta-data like sample alias, library, barcodes, etc...
      as well as read-level tags.) This tool takes an unaligned BAM with meta-data,
      and the aligned BAM produced by calling [SamToFastq](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_SamToFastq.php)
      and then passing the result to an aligner. It produces a new SAM file that includes
      all aligned and unaligned reads and also carries forward additional read attributes
      from the unmapped BAM (attributes that are otherwise lost in the process of
      converting to FASTQ). The resulting file will be valid for use by Picard and
      GATK tools. The output may be coordinate-sorted, in which case the tags, NM,
      MD, and UQ will be calculated and populated, or query-name sorted, in which
      case the tags will not be calculated or populated [1].\n\n*A list of **all inputs
      and parameters** with corresponding descriptions can be found at the bottom
      of the page.*\n\n###Common Use Cases\n\n* The **GATK MergeBamAlignment** tool
      requires a SAM or BAM file on its **Aligned BAM/SAM file** (`--ALIGNED_BAM`)
      input, original SAM or BAM file of unmapped reads, which must be in queryname
      order on its **Unmapped BAM/SAM file** (`--UNMAPPED_BAM`) input and a reference
      sequence on its **Reference** (`--REFERENCE_SEQUENCE`) input. The tool generates
      a single BAM/SAM file on its **Output merged BAM/SAM file** output.\n\n* Usage
      example:\n\n```\ngatk MergeBamAlignment \\\\\n      --ALIGNED_BAM aligned.bam
      \\\\\n      --UNMAPPED_BAM unmapped.bam \\\\\n      --OUTPUT merged.bam \\\\\n
      \     --REFERENCE_SEQUENCE reference_sequence.fasta\n```\n\n###Changes Introduced
      by Seven Bridges\n\n* The output file name will be prefixed using the **Output
      prefix** parameter. In case **Output prefix** is not provided, output prefix
      will be the same as the Sample ID metadata from **Input SAM/BAM file**, if the
      Sample ID metadata exists. Otherwise, output prefix will be inferred from the
      **Input SAM/BAM file** filename. This way, having identical names of the output
      files between runs is avoided. Moreover,  **merged** will be added before the
      extension of the output file name. \n\n* The user has a possibility to specify
      the output file format using the **Output file format** argument. Otherwise,
      the output file format will be the same as the format of the input aligned file.\n\n###Common
      Issues and Important Notes\n\n* Note:  This is not a tool for taking multiple
      BAM/SAM files and creating a bigger file by merging them. For that use-case,
      see [MergeSamFiles](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_MergeSamFiles.php).\n\n###Performance
      Benchmarking\n\nBelow is a table describing runtimes and task costs of **GATK
      MergeBamAlignment** for a couple of different samples, executed on the AWS cloud
      instances:\n\n| Experiment type |  Aligned BAM/SAM size |  Unmapped BAM/SAM
      size | Duration |  Cost | Instance (AWS) | \n|:--------------:|:------------:|:--------:|:-------:|:---------:|:----------:|:------:|:------:|------:|\n|
      \    RNA-Seq     |  1.4 GB |  1.9 GB |   9min   | ~0.06$ | c4.2xlarge (8 CPUs)
      | \n|     RNA-Seq     |  4.0 GB |  5.7 GB |   20min   | ~0.13$ | c4.2xlarge
      (8 CPUs) | \n|     RNA-Seq     | 6.6 GB | 9.5 GB |  32min  | ~0.21$ | c4.2xlarge
      (8 CPUs) | \n|     RNA-Seq     | 13 GB | 19 GB |  1h 4min  | ~0.42$ | c4.2xlarge
      (8 CPUs) |\n\n*Cost can be significantly reduced by using **spot instances**.
      Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances)
      for more details.*\n\n###References\n\n[1] [GATK MergeBamAlignment](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_MergeBamAlignment.php)"
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-mergebamalignment-4-1-0-0/14
    inputs:
    - doc: Adds the mate CIGAR tag (MC) if true, does not if false.
      id: add_mate_cigar
      inputBinding: {position: 4, prefix: --ADD_MATE_CIGAR, shellQuote: false}
      label: Add mate CIGAR
      sbg:altPrefix: -MC
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: add_mate_cigar
        symbols: ['true', 'false']
        type: enum
    - doc: Add PG tag to each read in a SAM or BAM.
      id: add_pg_tag_to_reads
      inputBinding: {position: 4, prefix: --ADD_PG_TAG_TO_READS, shellQuote: false}
      label: Add PG tag to reads
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: add_pg_tag_to_reads
        symbols: ['true', 'false']
        type: enum
    - doc: SAM or BAM file(s) with alignment data. Cannot be used in conjuction with
        argument(s) READ1_ALIGNED_BAM (R1_ALIGNED) READ2_ALIGNED_BAM (R2_ALIGNED).
      id: in_alignments
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   var arr = [].concat(inputs.in_alignments);\n    if (arr.length == 1)
          \n    {\n        return \"--ALIGNED_BAM \" + arr[0].path;\n    }\n    else\n
          \   {\n        var pe_1 = [];\n        var pe_2 = [];\n        var se =
          [];\n        for (var i in arr)\n        {\n            if (arr[i].metadata
          && arr[i].metadata.paired_end && arr[i].metadata.paired_end == 1)\n            {\n
          \               pe_1.push(arr[i].path);\n            }\n            else
          if (arr[i].metadata && arr[i].metadata.paired_end && arr[i].metadata.paired_end
          == 2)\n            {\n                pe_2.push(arr[i].path);\n            }\n
          \           else\n            {\n                se.push(arr[i].path);\n
          \           }\n        }\n        \n        if (se.length > 0) \n        {\n
          \           return \"--ALIGNED_BAM \" + se.join(\" --ALIGNED_BAM \");\n
          \       } \n        else if (pe_1.length > 0 && pe_2.length > 0 && pe_1.length
          == pe_2.length) \n        {\n            return \"--READ1_ALIGNED_BAM \"
          + pe_1.join(' --READ1_ALIGNED_BAM ') + \" --READ2_ALIGNED_BAM \" + pe_2.join('
          --READ2_ALIGNED_BAM ');\n        } \n        else \n        {\n            return
          \"\";\n        }\n            \n    }\n}"}
      label: Aligned BAM/SAM file
      sbg:category: Optional Arguments
      sbg:fileTypes: BAM, SAM
      sbg:toolDefaultValue: 'null'
      type: {items: File, type: array}
    - doc: Whether to output only aligned reads.
      id: aligned_reads_only
      inputBinding: {position: 4, prefix: --ALIGNED_READS_ONLY, shellQuote: false}
      label: Aligned reads only
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Use the aligner's idea of what a proper pair is rather than computing in
        this program.
      id: aligner_proper_pair_flags
      inputBinding: {position: 4, prefix: --ALIGNER_PROPER_PAIR_FLAGS, shellQuote: false}
      label: Aligner proper pair flags
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Attributes from the alignment record that should be removed when merging.
        This overrides ATTRIBUTES_TO_RETAIN if they share common tags.
      id: attributes_to_remove
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--ATTRIBUTES_TO_REMOVE',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Attributes to remove
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Reserved alignment attributes (tags starting with X, Y, or Z) that should
        be brought over from the alignment data when merging.
      id: attributes_to_retain
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--ATTRIBUTES_TO_RETAIN',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Attributes to retain
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Attributes on negative strand reads that need to be reversed.
      id: attributes_to_reverse
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--ATTRIBUTES_TO_REVERSE',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Attributes to reverse
      sbg:altPrefix: -RV
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '[OQ,U2]'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Attributes on negative strand reads that need to be reverse complemented.
      id: attributes_to_reverse_complement
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--ATTRIBUTES_TO_REVERSE_COMPLEMENT',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Attributes to reverse complement
      sbg:altPrefix: -RC
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '[E2,SQ]'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: Whether to clip adapters where identified.
      id: clip_adapters
      inputBinding: {position: 4, prefix: --CLIP_ADAPTERS, shellQuote: false}
      label: Clip adapters
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: clip_adapters
        symbols: ['true', 'false']
        type: enum
    - doc: For paired reads, soft clip the 3' end of each read if necessary so that
        it does not extend past the 5' end of its mate.
      id: clip_overlapping_reads
      inputBinding: {position: 4, prefix: --CLIP_OVERLAPPING_READS, shellQuote: false}
      label: Clip overlapping reads
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: clip_overlapping_reads
        symbols: ['true', 'false']
        type: enum
    - doc: Compression level for all compressed files created (e.g. BAM and VCF).
      id: compression_level
      inputBinding: {position: 4, prefix: --COMPRESSION_LEVEL, shellQuote: false}
      label: Compression level
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
      type: ['null', int]
    - doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
      id: create_index
      inputBinding: {position: 4, prefix: --CREATE_INDEX, shellQuote: false}
      label: Create index
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: The expected orientation of proper read pairs. Replaces JUMP_SIZE. Cannot
        be used in conjuction with argument(s) JUMP_SIZE (JUMP).
      id: expected_orientations
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--EXPECTED_ORIENTATIONS',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Expected orientations
      sbg:altPrefix: -ORIENTATIONS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: If false, do not write secondary alignments to output.
      id: include_secondary_alignments
      inputBinding: {position: 4, prefix: --INCLUDE_SECONDARY_ALIGNMENTS, shellQuote: false}
      label: Include secondary alignments
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: include_secondary_alignments
        symbols: ['true', 'false']
        type: enum
    - doc: Whether the lane is bisulfite sequence (used when calculating the NM tag).
      id: is_bisulfite_sequence
      inputBinding: {position: 4, prefix: --IS_BISULFITE_SEQUENCE, shellQuote: false}
      label: Is bisulfite sequence
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: The expected jump size (required if this is a jumping library). Deprecated.
        Use EXPECTED_ORIENTATIONS instead. Cannot be used in conjuction with argument(s)
        EXPECTED_ORIENTATIONS (ORIENTATIONS).
      id: jump_size
      inputBinding: {position: 4, prefix: --JUMP_SIZE, shellQuote: false}
      label: Jump size
      sbg:altPrefix: -JUMP
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', int]
    - doc: List of Sequence Records tags that must be equal (if present) in the reference
        dictionary and in the aligned file. Mismatching tags will cause an error if
        in this list, and a warning otherwise.
      id: matching_dictionary_tags
      inputBinding: {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n
          \   if (self)\n    {\n        var cmd = [];\n        for (var i = 0; i <
          self.length; i++) \n        {\n            cmd.push('--MATCHING_DICTIONARY_TAGS',
          self[i]);\n        }\n        return cmd.join(' ');\n    }\n    \n}"}
      label: Matching dictionary tags
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '[M5,LN]'
      type:
      - 'null'
      - {items: string, type: array}
    - doc: The maximum number of insertions or deletions permitted for an alignment
        to be included. Alignments with more than this many insertions or deletions
        will be ignored. Set to -1 to allow any number of insertions or deletions.
      id: max_insertions_or_deletions
      inputBinding: {position: 4, prefix: --MAX_INSERTIONS_OR_DELETIONS, shellQuote: false}
      label: Max insertions or deletions
      sbg:altPrefix: -MAX_GAPS
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    - doc: When writing files that need to be sorted, this will specify the number
        of records stored in RAM before spilling to disk. Increasing this number reduces
        the number of file handles needed to sort the file, and increases the amount
        of RAM needed.
      id: max_records_in_ram
      inputBinding: {position: 4, prefix: --MAX_RECORDS_IN_RAM, shellQuote: false}
      label: Max records in RAM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
      type: ['null', int]
    - doc: This input allows a user to set the desired overhead memory when running
        a tool or adding it to a workflow. This amount will be added to the Memory
        per job in the Memory requirements section but it will not be added to the
        -Xmx parameter leaving some memory not occupied which can be used as stack
        memory (-Xmx parameter defines heap memory). This input should be defined
        in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      id: memory_overhead_per_job
      label: Memory overhead per job
      sbg:category: Platform Options
      type: ['null', int]
    - doc: This input allows a user to set the desired memory requirement when running
        a tool or adding it to a workflow. This value should be propagated to the
        -Xmx parameter too.This input should be defined in MB (for both the platform
        part and the -Xmx part if Java tool is wrapped).
      id: memory_per_job
      label: Memory per job
      sbg:category: Platform Options
      type: ['null', int]
    - doc: If UNMAP_CONTAMINANT_READS is set, require this many unclipped bases or
        else the read will be marked as contaminant.
      id: min_unclipped_bases
      inputBinding: {position: 4, prefix: --MIN_UNCLIPPED_BASES, shellQuote: false}
      label: Min unclipped bases
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '32'
      type: ['null', int]
    - doc: DEPRECATED. This argument is ignored and will be removed.
      id: paired_run
      inputBinding: {position: 4, prefix: --PAIRED_RUN, shellQuote: false}
      label: Paired run
      sbg:altPrefix: -PE
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: paired_run
        symbols: ['true', 'false']
        type: enum
    - doc: 'Strategy for selecting primary alignment when the aligner has provided
        more than one alignment for a pair or fragment, and none are marked as primary,
        more than one is marked as primary, or the primary alignment is filtered out
        for some reason. For all strategies, ties are resolved arbitrarily. Possible
        values: { BestMapq (expects that multiple alignments will be correlated with
        HI tag, and prefers the pair of alignments with the largest MAPQ, in the absence
        of a primary selected by the aligner.) EarliestFragment (prefers the alignment
        which maps the earliest base in the read. Note that EarliestFragment may not
        be used for paired reads.) BestEndMapq (appropriate for cases in which the
        aligner is not pair-aware, and does not output the HI tag. It simply picks
        the alignment for each end with the highest MAPQ, and makes those alignments
        primary, regardless of whether the two alignments make sense together.) MostDistant
        (appropriate for a non-pair-aware aligner. Picks the alignment pair with the
        largest insert size. If all alignments would be chimeric, it picks the alignments
        for each end with the best MAPQ. ) }.'
      id: primary_alignment_strategy
      inputBinding: {position: 4, prefix: --PRIMARY_ALIGNMENT_STRATEGY, shellQuote: false}
      label: Primary alignment strategy
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: BestMapq
      type:
      - 'null'
      - name: primary_alignment_strategy
        symbols: [BestMapq, EarliestFragment, BestEndMapq, MostDistant]
        type: enum
    - doc: The command line of the program group (if not supplied by the aligned file).
      id: program_group_command_line
      inputBinding: {position: 4, prefix: --PROGRAM_GROUP_COMMAND_LINE, shellQuote: false}
      label: Program group command line
      sbg:altPrefix: -PG_COMMAND
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: The name of the program group (if not supplied by the aligned file).
      id: program_group_name
      inputBinding: {position: 4, prefix: --PROGRAM_GROUP_NAME, shellQuote: false}
      label: Program group name
      sbg:altPrefix: -PG_NAME
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: The version of the program group (if not supplied by the aligned file).
      id: program_group_version
      inputBinding: {position: 4, prefix: --PROGRAM_GROUP_VERSION, shellQuote: false}
      label: Program group version
      sbg:altPrefix: -PG_VERSION
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: The program group ID of the aligner (if not supplied by the aligned file).
      id: program_record_id
      inputBinding: {position: 4, prefix: --PROGRAM_RECORD_ID, shellQuote: false}
      label: Program record id
      sbg:altPrefix: -PG
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: The number of bases trimmed from the beginning of read 1 prior to alignment.
      id: read1_trim
      inputBinding: {position: 4, prefix: --READ1_TRIM, shellQuote: false}
      label: Read1 trim
      sbg:altPrefix: -R1_TRIM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: The number of bases trimmed from the beginning of read 2 prior to alignment.
      id: read2_trim
      inputBinding: {position: 4, prefix: --READ2_TRIM, shellQuote: false}
      label: Read2 trim
      sbg:altPrefix: -R2_TRIM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: Reference sequence file.
      id: in_reference
      inputBinding: {position: 4, prefix: --REFERENCE_SEQUENCE, shellQuote: false}
      label: Reference
      sbg:altPrefix: -R
      sbg:category: Required Arguments
      sbg:fileTypes: FASTA, FA
      secondaryFiles: [.fai, ^.dict]
      type: File
    - doc: The order in which the merged reads should be output.
      id: sort_order
      inputBinding: {position: 4, prefix: --SORT_ORDER, shellQuote: false}
      label: Sort order
      sbg:altPrefix: -SO
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: coordinate
      type:
      - 'null'
      - name: sort_order
        symbols: [unsorted, queryname, coordinate, duplicate, unknown]
        type: enum
    - doc: Detect reads originating from foreign organisms (e.g. bacterial DNA in
        a non-bacterial sample), and unmap + label those reads accordingly.
      id: unmap_contaminant_reads
      inputBinding: {position: 4, prefix: --UNMAP_CONTAMINANT_READS, shellQuote: false}
      label: Unmap contaminant reads
      sbg:altPrefix: -UNMAP_CONTAM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Original SAM or BAM file of unmapped reads, which must be in queryname
        order.
      id: unmapped_bam
      inputBinding: {position: 4, prefix: --UNMAPPED_BAM, shellQuote: false}
      label: Unmapped BAM/SAM file
      sbg:altPrefix: -UNMAPPED
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM
      type: File
    - doc: How to deal with alignment information in reads that are being unmapped
        (e.g. due to cross-species contamination.) Currently ignored unless UNMAP_CONTAMINANT_READS
        = true
      id: unmapped_read_strategy
      inputBinding: {position: 4, prefix: --UNMAPPED_READ_STRATEGY, shellQuote: false}
      label: Unmapped read strategy
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: DO_NOT_CHANGE
      type:
      - 'null'
      - name: unmapped_read_strategy
        symbols: [COPY_TO_TAG, DO_NOT_CHANGE, MOVE_TO_TAG]
        type: enum
    - doc: Validation stringency for all SAM files read by this program. Setting stringency
        to SILENT can improve performance when processing a BAM file in which variable-length
        data (read, qualities, tags) do not otherwise need to be decoded.
      id: validation_stringency
      inputBinding: {position: 4, prefix: --VALIDATION_STRINGENCY, shellQuote: false}
      label: Validation stringency
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
      type:
      - 'null'
      - name: validation_stringency
        symbols: [STRICT, LENIENT, SILENT]
        type: enum
    - doc: Output file name prefix.
      id: output_prefix
      label: Output prefix
      sbg:category: Optional Parameters
      type: ['null', string]
    - doc: Output file format
      id: output_file_format
      label: Output file format
      sbg:category: Optional parameters
      type:
      - 'null'
      - name: output_file_format
        symbols: [bam, sam]
        type: enum
    - doc: CPU per job.
      id: cpu_per_job
      label: CPU per job
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    label: GATK MergeBamAlignment
    outputs:
    - doc: Output merged SAM or BAM file.
      id: out_alignments
      label: Output merged SAM or BAM file
      outputBinding: {glob: '*am', outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: SAM, BAM
      secondaryFiles: ["${\n    if (self.nameext == \".bam\" && inputs.create_index)\n
          \   {\n        return [self.basename + \".bai\", self.nameroot + \".bai\"];\n
          \   }\n    else {\n        return []; \n    }\n}"]
      type: File
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: "${\n    return inputs.cpu_per_job ?
        inputs.cpu_per_job : 1;\n}", ramMin: "${\n    var memory = 4096;\n    if (inputs.memory_per_job)
        \n    {\n        memory = inputs.memory_per_job;\n    }\n    if (inputs.memory_overhead_per_job)\n
        \   {\n        memory += inputs.memory_overhead_per_job;\n    }\n    return
        memory;\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n
          \   else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar
          groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict
          = {};\n    for (var i = 0; i < files.length; i++) {\n        var value =
          files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};", "\nvar setMetadata = function(file,
          metadata) {\n    if (!('metadata' in file))\n        file['metadata'] =
          metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BAM Processing]
    sbg:content_hash: a758b43167e957642f45a0aad07716ff3b8c8c6a379cf76b35f10b0a3f5a121b8
    sbg:contributors: [uros_sipetic, nemanja.vucic, nens, veliborka_josipovic]
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/37
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552666475
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-mergebamalignment-4-1-0-0/14
    sbg:image_url: null
    sbg:latestRevision: 14
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - {id: 'https://software.broadinstitute.org/gatk/', label: Homepage}
    - {id: 'https://github.com/broadinstitute/gatk/', label: Source Code}
    - {id: 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip',
      label: Download}
    - {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
    - {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_MergeSamFiles.php',
      label: Documentation}
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1560336165
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 14
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/37
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1552666475, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/12}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554492767, 'sbg:revision': 1,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/23}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554720890, 'sbg:revision': 2,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/24}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554999266, 'sbg:revision': 3,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/25}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557734540, 'sbg:revision': 4, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/26}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000585, 'sbg:revision': 5, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/27}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558017849, 'sbg:revision': 6, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/28}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351570, 'sbg:revision': 7, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/29}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558370509, 'sbg:revision': 8, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/30}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558427482, 'sbg:revision': 9, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/31}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558448356, 'sbg:revision': 10, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/32}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558453788, 'sbg:revision': 11, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/33}
    - {'sbg:modifiedBy': nemanja.vucic, 'sbg:modifiedOn': 1559750464, 'sbg:revision': 12,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/34}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1560335266, 'sbg:revision': 13, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/36}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1560336165, 'sbg:revision': 14, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-mergebamalignment-4-1-0-0/37}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: -9
  sbg:y: 53.96965026855469
  scatter: [in_alignments, unmapped_bam]
  scatterMethod: dotproduct
- id: gatk_samtofastq_4_1_0_0
  in:
  - {default: true, id: include_non_pf_reads}
  - {id: in_alignments, source: in_alignments}
  - {default: true, id: interleave}
  label: GATK SamToFastq
  out:
  - {id: out_reads}
  - {id: unmapped_reads}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, prefix: '', shellQuote: false, valueFrom: "${\n    var in_alignments
        = [].concat(inputs.in_alignments)[0];\n    var output_ext    = inputs.compress_outputs
        ? \".fastq.gz\" : \".fastq\";\n    var interleave    = inputs.interleave;\n
        \   var output_prefix = ''; \n    var cmd_line      = '';\n    cmd_line          =
        \"cmd='' && paired_end=`samtools view -h \" + in_alignments.path + \" | head
        -n 500000 | samtools view -Sc -f 0x1 -`\";\n\n  if (!inputs.outputs_by_readgroup)\n
        \   {\n        if (inputs.output_prefix)\n        {\n            output_prefix
        = inputs.output_prefix;\n        }\n        else\n        {\n            if
        (in_alignments.metadata && in_alignments.metadata.sample_id)\n            {\n
        \               output_prefix = in_alignments.metadata.sample_id;\n            }\n
        \           else\n            {\n                output_prefix = in_alignments.path.split('/').pop().split('.')[0];\n
        \           }          \n        }           \n        \n        cmd_line
        = cmd_line + \" && if [ $paired_end != 0 ]; then cmd='--FASTQ \" + output_prefix;
        \n        \n        if (interleave)\n        {\n            cmd_line = cmd_line
        + \".interleaved\" + output_ext + \"';\";\n        }\n        else\n        {\n
        \           cmd_line = cmd_line + \".pe_1\" + output_ext;\n            cmd_line
        = cmd_line + \" --SECOND_END_FASTQ \" + output_prefix + \".pe_2\" + output_ext;\n
        \           cmd_line = cmd_line + \" --UNPAIRED_FASTQ \" + output_prefix +
        \".unpaired\" + output_ext + \"';\";\n        }        \n        cmd_line
        = cmd_line + \" else cmd='--FASTQ \" + output_prefix  + \".se\" + output_ext
        + \"'; fi;\";\n        return cmd_line;\n    }\n    else\n    {\n        return
        \"cmd='--OUTPUT_DIR .'\";\n    }\n}\n\n"}
    - {position: 1, prefix: '', shellQuote: false, valueFrom: /opt/gatk}
    - {position: 2, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)\n
        \   {\n        return \"--java-options\";\n    }\n    else {\n        return
        '';\n    }\n    \n}"}
    - {position: 3, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)
        {\n        return '\\\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\\\"';\n
        \   }\n    else {\n        return \"\";\n    }\n}"}
    - {position: 4, prefix: '', shellQuote: false, valueFrom: SamToFastq}
    - {position: 6, prefix: '', shellQuote: false, valueFrom: "${\n        return
        '$cmd';\n}"}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "The **GATK SamToFastq** tool converts a SAM or BAM file to FASTQ.\n\nThis
      tool extracts read sequences and qualities from the input SAM/BAM file and writes
      them into the output file in Sanger FASTQ format.\n\nIn the RC mode (default
      is True), if the read is aligned and the alignment is to the reverse strand
      on the genome, the read sequence from input SAM file will be reverse-complemented
      prior to writing it to FASTQ in order to correctly restore the original read
      sequence as it was generated by the sequencer [1].\n\n*A list of **all inputs
      and parameters** with corresponding descriptions can be found at the bottom
      of the page.*\n\n###Common Use Cases\n\n* The **GATK SamToFastq** tool requires
      a BAM/SAM file on its **Input BAM/SAM file** (`--INPUT`) input. The tool generates
      a single-end FASTQ file on its **Output FASTQ file(s)** output if the input
      BAM/SAM file is single end. In case the input file is paired end, the tool outputs
      the first end of the pair FASTQ and the second end of the pair FASTQ on its
      **Output FASTQ file(s)** output, except when the **Interleave** (`--INTERLEAVE`)
      option is set to True. If the output is an interleaved FASTQ file, if paired,
      each line will have /1 or /2 to describe which end it came from.\n\n* The **GATK
      SamToFastq** tool supports an optional parameter  **Output by readgroup** (`--OUTPUT_BY_READGROUP`)
      which, when true, outputs a FASTQ file per read group (two FASTQ files per read
      group if the group is paired).\n\n* Usage example (input BAM file is single-end):\n\n```\ngatk
      SamToFastq \n     --INPUT input.bam\n     --FASTQ output.fastq\n```\n\n\n\n\n\n*
      Usage example (input BAM file is paired-end):\n\n```\ngatk SamToFastq \n     --INPUT
      input.bam\n     --FASTQ output.pe_1.fastq\n     --SECOND_END_FASTQ output.pe_2.fastq\n
      \    --UNPAIRED_FASTQ unpaired.fastq\n\n```\n\n###Changes Introduced by Seven
      Bridges\n\n* The GATK SamToFastq tool is implemented to check if the input alignments
      file contains single-end or paired-end data and according to that generates
      different command lines for these two modes and thus produces appropriate output
      files on its **Output FASTQ file(s)** output (one FASTQ file in single-end mode
      and two FASTQ files if the input alignment file contains paired-end data). \n\n*
      All output files will be prefixed using the **Output prefix** parameter. In
      case the **Output prefix** is not provided, the output prefix will be the same
      as the Sample ID metadata from the **input SAM/BAM file**, if the Sample ID
      metadata exists. Otherwise, the output prefix will be inferred from the **Input
      SAM/BAM** filename. This way, having identical names of the output files between
      runs is avoided.\n\n* For paired-end read files, in order to successfully run
      alignment with STAR, this tool adds the appropriate **paired-end** metadata
      field in the output FASTQ files.\n\n###Common Issues and Important Notes\n\n*
      None\n\n###Performance Benchmarking\n\nBelow is a table describing runtimes
      and task costs of **GATK SamToFastq** for a couple of different samples, executed
      on the AWS cloud instances:\n\n| Experiment type |  Input size | Paired-end
      | # of reads | Read length | Duration |  Cost | Instance (AWS) | \n|:--------------:|:------------:|:--------:|:-------:|:---------:|:----------:|:------:|:------:|\n|
      \    RNA-Seq     |  1.9 GB |     Yes    |     16M     |     101     |   4min
      \  | ~0.03$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq     |  5.7 GB |     Yes
      \   |     50M     |     101     |   7min   | ~0.04$ | c4.2xlarge (8 CPUs) |
      \n|     RNA-Seq     | 9.5 GB |     Yes    |     82M    |     101     |  10min
      \ | ~0.07$ | c4.2xlarge (8 CPUs) | \n|     RNA-Seq     | 19 GB |     Yes    |
      \    164M    |     101     |  20min  | ~0.13$ | c4.2xlarge (8 CPUs) |\n\n*Cost
      can be significantly reduced by using **spot instances**. Visit the [Knowledge
      Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n\n###References\n\n[1]
      [GATK SamToFastq](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_sam_SamToFastq)"
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-samtofastq-4-1-0-0/15
    inputs:
    - doc: 'The action that should be taken with clipped reads: ''X'' means the reads
        and qualities should be trimmed at the clipped position; ''N'' means the bases
        should be changed to Ns in the clipped region; and any integer means that
        the base qualities should be set to that value in the clipped region.'
      id: clipping_action
      inputBinding: {position: 5, prefix: --CLIPPING_ACTION, shellQuote: false}
      label: Clipping action
      sbg:altPrefix: -CLIP_ACT
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: The attribute that stores the position at which the SAM record should be
        clipped.
      id: clipping_attribute
      inputBinding: {position: 5, prefix: --CLIPPING_ATTRIBUTE, shellQuote: false}
      label: Clipping attribute
      sbg:altPrefix: -CLIP_ATTR
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', string]
    - doc: When performing clipping with the CLIPPING_ATTRIBUTE and CLIPPING_ACTION
        parameters, ensure that the resulting reads after clipping are at least CLIPPING_MIN_LENGTH
        bases long. If the original read is shorter than CLIPPING_MIN_LENGTH then
        the original read length will be maintained.
      id: clipping_min_length
      inputBinding: {position: 5, prefix: --CLIPPING_MIN_LENGTH, shellQuote: false}
      label: Clipping min length
      sbg:altPrefix: -CLIP_MIN
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: Compress output FASTQ files per read group using gzip and append a .gz
        extension to the file names. Cannot be used in conjuction with argument(s)
        FASTQ (F) SECOND_END_FASTQ (F2) UNPAIRED_FASTQ (FU).
      id: compress_outputs_per_rg
      inputBinding: {position: 5, prefix: --COMPRESS_OUTPUTS_PER_RG, shellQuote: false}
      label: Compress outputs per RG
      sbg:altPrefix: -GZOPRG
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Compression level for all compressed files created (e.g. BAM and VCF).
      id: compression_level
      inputBinding: {position: 5, prefix: --COMPRESSION_LEVEL, shellQuote: false}
      label: Compression level
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
      type: ['null', int]
    - doc: Include non-PF reads from the SAM file into the output FASTQ files. PF
        means 'passes filtering'. Reads whose 'not passing quality controls' flag
        is set are non-PF reads. See GATK Dictionary for more info.
      id: include_non_pf_reads
      inputBinding: {position: 5, prefix: --INCLUDE_NON_PF_READS, shellQuote: false}
      label: Include non PF reads
      sbg:altPrefix: -NON_PF
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: If true, include non-primary alignments in the output. Support of non-primary
        alignments in SamToFastq is not comprehensive, so there may be exceptions
        if this is set to true and there are paired reads with non-primary alignments.
      id: include_non_primary_alignments
      inputBinding: {position: 5, prefix: --INCLUDE_NON_PRIMARY_ALIGNMENTS, shellQuote: false}
      label: Include non primary alignments
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Input SAM/BAM file to extract reads from.
      id: in_alignments
      inputBinding: {position: 5, prefix: --INPUT, shellQuote: false}
      label: Input SAM/BAM file
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: SAM, BAM
      type: File
    - doc: Will generate an interleaved FASTQ if paired, each line will have /1 or
        /2 to describe which end it came from.
      id: interleave
      inputBinding: {position: 5, prefix: --INTERLEAVE, shellQuote: false}
      label: Interleave
      sbg:altPrefix: -INTER
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: This input allows a user to set the desired overhead memory when running
        a tool or adding it to a workflow. This amount will be added to the Memory
        per job in the Memory requirements section but it will not be added to the
        -Xmx parameter leaving some memory not occupied which can be used as stack
        memory (-Xmx parameter defines heap memory). This input should be defined
        in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      id: memory_overhead_per_job
      label: Memory overhead per job
      sbg:category: Platform Options
      type: ['null', int]
    - doc: This input allows a user to set the desired memory requirement when running
        a tool or adding it to a workflow. This value should be propagated to the
        -Xmx parameter too.This input should be defined in MB (for both the platform
        part and the -Xmx part if Java tool is wrapped).
      id: memory_per_job
      label: Memory per job
      sbg:category: Platform Options
      sbg:toolDefaultValue: 2048 MB
      type: ['null', int]
    - doc: Output a FASTQ file per read group (two FASTQ files per read group if the
        group is paired). Cannot be used in conjuction with argument(s)FASTQ (F) SECOND_END_FASTQ
        (F2) UNPAIRED_FASTQ (FU).
      id: output_per_rg
      inputBinding: {position: 5, prefix: --OUTPUT_PER_RG, shellQuote: false}
      label: Output per RG
      sbg:altPrefix: -OPRG
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: End-trim reads using the phred/bwa quality trimming algorithm and this
        quality.
      id: quality
      inputBinding: {position: 5, prefix: --QUALITY, shellQuote: false}
      label: Quality
      sbg:altPrefix: -Q
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', int]
    - doc: Re-reverse bases and qualities of reads with negative strand flag set before
        writing them to FASTQ.
      id: re_reverse
      inputBinding: {position: 5, prefix: --RE_REVERSE, shellQuote: false}
      label: Re reverse
      sbg:altPrefix: -RC
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'true'
      type:
      - 'null'
      - name: re_reverse
        symbols: ['true', 'false']
        type: enum
    - doc: The maximum number of bases to write from read 1 after trimming. If there
        are fewer than this many bases left after trimming, all will be written. If
        this value is null then all bases left after trimming will be written.
      id: read1_max_bases_to_write
      inputBinding: {position: 5, prefix: --READ1_MAX_BASES_TO_WRITE, shellQuote: false}
      label: Read1 max bases to write
      sbg:altPrefix: -R1_MAX_BASES
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', int]
    - doc: The number of bases to trim from the beginning of read 1.
      id: read1_trim
      inputBinding: {position: 5, prefix: --READ1_TRIM, shellQuote: false}
      label: Read1 trim
      sbg:altPrefix: -R1_TRIM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: The maximum number of bases to write from read 2 after trimming. If there
        are fewer than this many bases left after trimming, all will be written. If
        this value is null then all bases left after trimming will be written.
      id: read2_max_bases_to_write
      inputBinding: {position: 5, prefix: --READ2_MAX_BASES_TO_WRITE, shellQuote: false}
      label: Read2 max bases to write
      sbg:altPrefix: -R2_MAX_BASES
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'null'
      type: ['null', int]
    - doc: The number of bases to trim from the beginning of read 2.
      id: read2_trim
      inputBinding: {position: 5, prefix: --READ2_TRIM, shellQuote: false}
      label: Read2 trim
      sbg:altPrefix: -R2_TRIM
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '0'
      type: ['null', int]
    - doc: The read group tag (PU or ID) to be used to output a FASTQ file per read
        group.
      id: rg_tag
      inputBinding: {position: 5, prefix: --RG_TAG, shellQuote: false}
      label: RG tag
      sbg:altPrefix: -RGT
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: PU
      type: ['null', string]
    - doc: Validation stringency for all SAM files read by this program. Setting stringency
        to silent can improve performance when processing a BAM file in which variable-length
        data (read, qualities, tags) do not otherwise need to be decoded.
      id: validation_stringency
      inputBinding: {position: 5, prefix: --VALIDATION_STRINGENCY, shellQuote: false}
      label: Validation stringency
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
      type:
      - 'null'
      - name: validation_stringency
        symbols: [STRICT, LENIENT, SILENT]
        type: enum
    - doc: Output file name prefix.
      id: output_prefix
      label: Output prefix
      sbg:category: Optional Arguments
      type: ['null', string]
    - doc: Compress output file(s).
      id: compress_outputs
      label: Compress output file(s)
      sbg:category: Optional parameters
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: CPU per job.
      id: cpu_per_job
      label: CPU per job
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    label: GATK SamToFastq
    outputs:
    - doc: Output FASTQ file(s).
      id: out_reads
      label: Output FASTQ file(s)
      outputBinding: {glob: "${\n    var output_ext = inputs.compress_outputs ? \".fastq.gz\"
          : \".fastq\";\n    var interleave = inputs.interleave;\n    if (!inputs.outputs_by_readgroup)\n
          \   {\n        if (interleave)\n        {\n            return \"*interleaved\"
          + output_ext;\n        }\n        else\n        {\n            return [\"*pe_1\"
          + output_ext, \"*pe_2\" + output_ext, \"*se\" + output_ext];\n        }\n\n
          \   }\n    else\n    {\n        return \"*\" + output_ext;\n    }\n}", outputEval: "${
          \n    self = [].concat(self)\n    \n    function getPairedEnd(filename)\n
          \   {\n        if (filename.lastIndexOf(\".fastq\") !== 0 && filename[filename.lastIndexOf(\".fastq\")
          - 2 ]==\"_\") \n        {\n            return filename[filename.lastIndexOf(\".fastq\")
          - 1 ];\n        } \n        else \n        {\n            return \"\";\n
          \       }\n    }\n    \n    var out = inheritMetadata(self,inputs.in_alignments);\n
          \   for (var i=0; i < out.length; i++)\n    {\n        out[i].metadata['paired_end']
          = getPairedEnd(out[i].path);\n    }\n    \n    return out;\n}"}
      sbg:fileTypes: FASTQ, FASTQ.GZ
      type:
      - 'null'
      - {items: File, type: array}
    - doc: Unpaired reads.
      id: unmapped_reads
      label: Unpaired reads
      outputBinding: {glob: "${\n    var output_ext = inputs.compress_outputs ? \".fastq.gz\"
          : \".fastq\";\n    var interleave = inputs.interleave;\n    if (!inputs.outputs_by_readgroup)\n
          \   {\n        if (!interleave)\n        {\n            return \"*unpaired\"
          + output_ext;\n        }\n        else \n        {\n             return
          \"\"; \n        }      \n    }\n  else {\n       return \"\";  \n  \n  }\n}",
        outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: FASTQ, FASTQ.GZ
      type: ['null', File]
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: "${\n    return inputs.cpu_per_job ?
        inputs.cpu_per_job : 1;\n}", ramMin: "${\n    var memory = 4096;\n    if (inputs.memory_per_job)
        \n    {\n        memory = inputs.memory_per_job;\n    }\n    if (inputs.memory_overhead_per_job)\n
        \   {\n        memory += inputs.memory_overhead_per_job;\n    }\n    return
        memory;\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n
          \   else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar
          groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict
          = {};\n    for (var i = 0; i < files.length; i++) {\n        var value =
          files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};", "\nvar setMetadata = function(file,
          metadata) {\n    if (!('metadata' in file)) {\n        file['metadata']
          = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key]
          = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1,
          o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2
          = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example
          = o2[i]['metadata'];\n        for (var key in example) {\n            if
          (i == 0)\n                commonMetadata[key] = example[key];\n            else
          {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete
          commonMetadata[key]\n                }\n            }\n        }\n    }\n
          \   if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n
          \   } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i]
          = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BAM Processing]
    sbg:content_hash: a21e1194e724a1f17bceabd4d2040324713c2a5c63896adcebbc777578b2bfef5
    sbg:contributors: [uros_sipetic, nens, veliborka_josipovic]
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1552663400
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-samtofastq-4-1-0-0/15
    sbg:image_url: null
    sbg:latestRevision: 15
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - {id: 'https://software.broadinstitute.org/gatk/', label: Homepage}
    - {id: 'https://github.com/broadinstitute/gatk/', label: Source Code}
    - {id: 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip',
      label: Download}
    - {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
    - {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_SamToFastq.php',
      label: Documentation}
    sbg:modifiedBy: veliborka_josipovic
    sbg:modifiedOn: 1561548030
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 15
    sbg:revisionNotes: Added glob for single end output fastq
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1552663400, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/19}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1552663734, 'sbg:revision': 1,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/20}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554492676, 'sbg:revision': 2,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/27}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554493243, 'sbg:revision': 3,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/28}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554720826, 'sbg:revision': 4,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/29}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1554999298, 'sbg:revision': 5,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-samtofastq-4-1-0-0/30}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557484228, 'sbg:revision': 6, 'sbg:revisionNotes': Updated
        Description}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557745933, 'sbg:revision': 7, 'sbg:revisionNotes': ''}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557918579, 'sbg:revision': 8, 'sbg:revisionNotes': 'v32:
        [input]'}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557919927, 'sbg:revision': 9, 'sbg:revisionNotes': v5->update}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558002424, 'sbg:revision': 10, 'sbg:revisionNotes': output
        required}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558015833, 'sbg:revision': 11, 'sbg:revisionNotes': unmapped_reads
        required}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558022954, 'sbg:revision': 12, 'sbg:revisionNotes': strict
        js for glob}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558023482, 'sbg:revision': 13, 'sbg:revisionNotes': strict
        javascript for unmapped_reads}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558354170, 'sbg:revision': 14, 'sbg:revisionNotes': return
        '';}
    - {'sbg:modifiedBy': veliborka_josipovic, 'sbg:modifiedOn': 1561548030, 'sbg:revision': 15,
      'sbg:revisionNotes': Added glob for single end output fastq}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: -440.40716552734375
  sbg:y: 290.9440612792969
  scatter: [in_alignments]
- id: gatk_setnmmdanduqtags_4_1_0_0
  in:
  - {default: true, id: create_index}
  - {id: in_alignments, source: gatk_sortsam_4_1_0_0/out_alignments}
  - {id: reference_sequence, source: in_reference}
  label: GATK SetNmMdAndUqTags
  out:
  - {id: out_alignments}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, shellQuote: false, valueFrom: /opt/gatk}
    - {position: 1, shellQuote: false, valueFrom: --java-options}
    - {position: 2, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)
        {\n        return '\\\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\\\"';\n
        \   }\n    return '\\\"-Xmx2048M\\\"';\n}"}
    - {position: 3, shellQuote: false, valueFrom: SetNmMdAndUqTags}
    - {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n    var tmp = [].concat(inputs.in_alignments);\n
        \   var ext = \"\"; \n    if (inputs.output_file_format) {\n        ext =
        inputs.output_file_format;\n    } else {\n        ext = tmp[0].path.split('.').pop();\n
        \   }\n    \n    if (inputs.output_prefix) {\n        return '-O ' +  inputs.output_prefix
        + \".fixed.\" + ext;\n    } else if (tmp[0].metadata && tmp[0].metadata.sample_id)
        {\n        return '-O ' +  tmp[0].metadata.sample_id + \".fixed.\" + ext;\n
        \   } else {\n        return '-O ' +  tmp[0].path.split('/').pop().split(\".\")[0]
        + \".fixed.\" + ext;\n    }\n    \n}"}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "The **GATK SetNmMdAndUqTags** tool takes in a coordinate-sorted SAM or BAM
      and calculatesthe NM, MD, and UQ tags by comparing it with the reference. \n\nThe
      **GATK SetNmMdAndUqTags**  may be needed when **GATK MergeBamAlignment** was
      run with **SORT_ORDER** other than `coordinate` and thus could not fix these
      tags. \n\n\n###Common Use Cases\nThe **GATK SetNmMdAndUqTags** tool  fixes NM,
      MD and UQ tags in SAM/BAM file **Input SAM/BAM file**   (`--INPUT`)  input.
      This tool takes in a coordinate-sorted SAM or BAM file and calculates the NM,
      MD, and UQ tags by comparing with the reference **Reference sequence** (`--REFERENCE_SEQUENCE`).\n\n*
      Usage example:\n\n```\njava -jar picard.jar SetNmMdAndUqTags\n     --REFERENCE_SEQUENCE=reference_sequence.fasta\n
      \    --INPUT=sorted.bam\n```\n\n\n###Changes Introduced by Seven Bridges\n\n*
      Prefix of the output file is defined with the optional parameter **Output prefix**.
      If **Output prefix** is not provided, name of the sorted file is obtained from
      **Sample ID** metadata form the **Input SAM/BAM file**, if the **Sample ID**
      metadata exists. Otherwise, the output prefix will be inferred form the **Input
      SAM/BAM file** filename. \n\n\n\n###Common Issues and Important Notes\n\n* The
      **Input SAM/BAM file** must be coordinate sorted in order to run  **GATK SetNmMdAndUqTags**.
      \n* If specified, the MD and NM tags can be ignored and only the UQ tag be set.
      \n\n\n###References\n[1] [GATK SetNmMdAndUqTags home page](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_sam_SetNmMdAndUqTags.php)"
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-setnmmdanduqtags-4-1-0-0/10
    inputs:
    - doc: Validation stringency for all sam files read by this program. Setting stringency
        to silent can improve performance when processing a bam file in which variable-length
        data (read, qualities, tags) do not otherwise need to be decoded.
      id: validation_stringency
      inputBinding: {position: 4, prefix: --VALIDATION_STRINGENCY, shellQuote: false}
      label: Validation stringency
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
      type:
      - 'null'
      - name: validation_stringency
        symbols: [STRICT, LENIENT, SILENT]
        type: enum
    - doc: The fixed bam or sam output prefix name.
      id: output_prefix
      label: Output
      sbg:altPrefix: -O
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: sample_id.fixed.bam
      type: ['null', string]
    - doc: This input allows a user to set the desired overhead memory when running
        a tool or adding it to a workflow. This amount will be added to the Memory
        per job in the Memory requirements section but it will not be added to the
        -Xmx parameter leaving some memory not occupied which can be used as stack
        memory (-Xmx parameter defines heap memory). This input should be defined
        in MB (for both the platform part and the -Xmx part if Java tool is wrapped).
      id: memory_overhead_per_job
      label: Memory Overhead Per Job
      sbg:category: Execution
      type: ['null', int]
    - doc: When writing files that need to be sorted, this will specify the number
        of records stored in ram before spilling to disk. Increasing this number reduces
        the number of file handles needed to sort the file, and increases the amount
        of ram needed.
      id: max_records_in_ram
      inputBinding: {position: 4, prefix: --MAX_RECORDS_IN_RAM, shellQuote: false}
      label: Max records in ram
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
      type: ['null', int]
    - doc: Whether to create a bam index when writing a coordinate-sorted bam file.
      id: create_index
      inputBinding: {position: 4, prefix: --CREATE_INDEX, shellQuote: false}
      label: Create index
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Whether the file contains bisulfite sequence (used when calculating the
        nm tag).
      id: is_bisulfite_sequence
      inputBinding: {position: 4, prefix: --IS_BISULFITE_SEQUENCE, shellQuote: false}
      label: Is bisulfite sequence
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Compression level for all compressed files created (e.g. Bam and vcf).
      id: compression_level
      inputBinding: {position: 4, prefix: --COMPRESSION_LEVEL, shellQuote: false}
      label: Compression level
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
      type: ['null', int]
    - doc: This input allows a user to set the desired memory requirement when running
        a tool or adding it to a workflow. This value should be propagated to the
        -Xmx parameter too.This input should be defined in MB (for both the platform
        part and the -Xmx part if Java tool is wrapped).
      id: memory_per_job
      label: Memory Per Job
      sbg:category: Execution
      type: ['null', int]
    - doc: The BAM or SAM file to fix.
      id: in_alignments
      inputBinding: {position: 4, prefix: --INPUT, shellQuote: false}
      label: Input SAM/BAM file
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM
      type: File
    - doc: Reference sequence FASTA file.
      id: reference_sequence
      inputBinding: {position: 4, prefix: --REFERENCE_SEQUENCE, shellQuote: false}
      label: Reference sequence
      sbg:altPrefix: -R
      sbg:category: Required Arguments
      sbg:fileTypes: FASTA, FA
      type: File
    - doc: Only set the uq tag, ignore md and nm.
      id: set_only_uq
      inputBinding: {position: 4, prefix: --SET_ONLY_UQ, shellQuote: false}
      label: Set only uq
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: This input allows a user to set the desired CPU requirement when running
        a tool or adding it to a workflow.
      id: cpu_per_job
      label: CPU per job
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    - doc: Output file format.
      id: output_file_format
      label: Output file format
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: Same as input
      type:
      - 'null'
      - name: output_file_format
        symbols: [bam, sam]
        type: enum
    label: GATK SetNmMdAndUqTags
    outputs:
    - doc: Output BAM/SAM file with fixed tags.
      id: out_alignments
      label: Output BAM/SAM file
      outputBinding: {glob: '*am', outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: BAM, SAM
      secondaryFiles: ["${  \n    if (inputs.create_index)\n    {\n        return
          self.nameroot + \".bai\";\n    }\n    else {\n        return ''; \n    }\n}"]
      type: ['null', File]
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: "${\n    return inputs.cpu_per_job ?
        inputs.cpu_per_job : 1;\n}", ramMin: "${\n    var memory = 4096;\n    if (inputs.memory_per_job)
        \n    {\n        memory = inputs.memory_per_job;\n    }\n    if (inputs.memory_overhead_per_job)\n
        \   {\n        memory += inputs.memory_overhead_per_job;\n    }\n    return
        memory;\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n
          \   for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n
          \   }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n
          \   var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2
          = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example
          = o2[i]['metadata'];\n        for (var key in example) {\n            if
          (i == 0)\n                commonMetadata[key] = example[key];\n            else
          {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete
          commonMetadata[key]\n                }\n            }\n        }\n    }\n
          \   if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n
          \   } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i]
          = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar
          toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy
          = function(files, key) {\n    var groupedFiles = [];\n    var tempDict =
          {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n
          \       if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};", "\nvar setMetadata = function(file,
          metadata) {\n    if (!('metadata' in file))\n        file['metadata'] =
          metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BAM Processing]
    sbg:content_hash: a31d48359c8ea5e8ac91b2096488ac9e8a71d49dd3aa1a8ffbdcc09665a2c1f39
    sbg:contributors: [nens, uros_sipetic]
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/15
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1555498307
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-setnmmdanduqtags-4-1-0-0/10
    sbg:image_url: null
    sbg:latestRevision: 10
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - {id: 'https://software.broadinstitute.org/gatk/', label: Homepage}
    - {id: 'https://github.com/broadinstitute/gatk/', label: Source Code}
    - {id: 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip',
      label: Download}
    - {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
    - {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_SetNmMdAndUqTags.php',
      label: Documentation}
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1558518048
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 10
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/15
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1555498307, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/1}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1555582274, 'sbg:revision': 1,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/5}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1556194603, 'sbg:revision': 2,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/6}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557399646, 'sbg:revision': 3, 'sbg:revisionNotes': app
        info improved - perf bench needed}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557417063, 'sbg:revision': 4, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/7}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557734531, 'sbg:revision': 5, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/9}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000576, 'sbg:revision': 6, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/10}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558100350, 'sbg:revision': 7, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/11}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351574, 'sbg:revision': 8, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/13}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558450064, 'sbg:revision': 9, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/14}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558518048, 'sbg:revision': 10, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-setnmmdanduqtags-4-1-0-0/15}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: 617.387451171875
  sbg:y: 223.81729125976562
- id: gatk_sortsam_4_1_0_0
  in:
  - {id: in_alignments, source: gatk_markduplicates_4_1_0_0/out_alignments}
  - {default: coordinate, id: sort_order}
  label: GATK SortSam
  out:
  - {id: out_alignments}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    arguments:
    - {position: 0, shellQuote: false, valueFrom: /opt/gatk}
    - {position: 1, shellQuote: false, valueFrom: --java-options}
    - {position: 2, prefix: '', shellQuote: false, valueFrom: "${\n    if (inputs.memory_per_job)
        {\n        return '\\\"-Xmx'.concat(inputs.memory_per_job, 'M') + '\\\"';\n
        \   }\n    return '\\\"-Xmx2048M\\\"';\n}"}
    - {position: 3, shellQuote: false, valueFrom: SortSam}
    - {position: 4, prefix: '', shellQuote: false, valueFrom: "${\n    var tmp = [].concat(inputs.in_alignments);\n
        \   var ext = '';\n  \n    if (inputs.output_file_format){\n        ext =
        inputs.output_file_format;\n    }    else {\n        ext = tmp[0].path.split(\".\").pop();\n
        \   }\n    \n    \n    if (inputs.output_prefix) {\n        return '-O ' +
        \ inputs.output_prefix + \".sorted.\" + ext;\n      \n    }else if (tmp[0].metadata
        && tmp[0].metadata.sample_id) {\n        \n        return '-O ' +  tmp[0].metadata.sample_id
        + \".sorted.\" + ext;\n    } else {\n         \n        return '-O ' +  tmp[0].path.split('/').pop().split(\".\")[0]
        + \".sorted.\"+ext;\n    }\n    \n    \n}"}
    baseCommand: []
    class: CommandLineTool
    cwlVersion: v1.0
    doc: "The **GATK SortSam** tool sorts the input SAM or BAM file by coordinate,
      queryname (QNAME), or some other property of the SAM record.\n\nThe **GATK SortOrder**
      of a SAM/BAM file is found in the SAM file header tag @HD in the field labeled
      SO.  For a coordinate\nsorted SAM/BAM file, read alignments are sorted first
      by the reference sequence name (RNAME) field using the reference\nsequence dictionary
      (@SQ tag).  Alignments within these subgroups are secondarily sorted using the
      left-most mapping\nposition of the read (POS).  Subsequent to this sorting scheme,
      alignments are listed arbitrarily.</p><p>For\nqueryname-sorted alignments, all
      alignments are grouped using the queryname field but the alignments are not
      necessarily\nsorted within these groups.  Reads having the same queryname are
      derived from the same template\n\n\n###Common Use Cases\n\nThe **GATK SortSam**
      tool requires a BAM/SAM file on its **Input SAM/BAM file**   (`--INPUT`)  input.
      The tool sorts input file in the order defined by (`--SORT_ORDER`) parameter.
      Available sort order options are `queryname`, `coordinate` and `duplicate`.
      \ \n\n* Usage example:\n\n```\njava -jar picard.jar SortSam\n     --INPUT=input.bam
      \n     --SORT_ORDER=coordinate\n```\n\n\n###Changes Introduced by Seven Bridges\n\n*
      Prefix of the output file is defined with the optional parameter **Output prefix**.
      If **Output prefix** is not provided, name of the sorted file is obtained from
      **Sample ID** metadata from the **Input SAM/BAM file**, if the **Sample ID**
      metadata exists. Otherwise, the output prefix will be inferred form the **Input
      SAM/BAM file** filename. \n\n\n###Common Issues and Important Notes\n\n* None\n\n\n###Performance
      Benchmarking\nBelow is a table describing runtimes and task costs of **GATK
      SortSam** for a couple of different samples, executed on the AWS cloud instances:\n\n|
      Experiment type |  Input size | Paired-end | # of reads | Read length | Duration
      |  Cost | Instance (AWS) | \n|:--------------:|:------------:|:--------:|:-------:|:---------:|:----------:|:------:|:------:|\n|
      \    WGS     |          |     Yes    |     16M     |     101     |   4min   |
      ~0.03$ | c4.2xlarge (8 CPUs) | \n|     WGS     |         |     Yes    |     50M
      \    |     101     |   7min   | ~0.04$ | c4.2xlarge (8 CPUs) | \n|     WGS     |
      \        |     Yes    |     82M    |     101     |  10min  | ~0.07$ | c4.2xlarge
      (8 CPUs) | \n|     WES     |         |     Yes    |     164M    |     101     |
      \ 20min  | ~0.13$ | c4.2xlarge (8 CPUs) |\n\n*Cost can be significantly reduced
      by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances)
      for more details.*\n\n\n\n###References\n[1] [GATK SortSam home page](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_sam_SortSam.php)"
    id: uros_sipetic/gatk-4-1-0-0-demo/gatk-sortsam-4-1-0-0/8
    inputs:
    - doc: Input BAM or SAM file to sort.  Required
      id: in_alignments
      inputBinding: {position: 4, prefix: --INPUT, shellQuote: false}
      label: Input SAM/BAM file
      sbg:altPrefix: -I
      sbg:category: Required Arguments
      sbg:fileTypes: BAM, SAM
      type: File
    - doc: Sorted bam or sam output file.
      id: output_prefix
      label: Output prefix
      sbg:altPrefix: -O
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: sample_id.sorted.bam
      type: ['null', string]
    - doc: Compression level for all compressed files created (e.g. Bam and vcf).
      id: compression_level
      inputBinding: {position: 4, prefix: --COMPRESSION_LEVEL, shellQuote: false}
      label: Compression level
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '2'
      type: ['null', int]
    - doc: Whether to create a bam index when writing a coordinate-sorted bam file.
      id: create_index
      inputBinding: {position: 4, prefix: --CREATE_INDEX, shellQuote: false}
      label: Create index
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: Whether to create an md5 digest for any bam or fastq files created.
      id: create_md5_file
      inputBinding: {position: 4, prefix: --CREATE_MD5_FILE, shellQuote: false}
      label: Create md5 file
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: 'false'
      type: ['null', boolean]
    - doc: When writing files that need to be sorted, this will specify the number
        of records stored in ram before spilling to disk. Increasing this number reduces
        the number of file handles needed to sort the file, and increases the amount
        of ram needed.
      id: max_records_in_ram
      inputBinding: {position: 4, prefix: --MAX_RECORDS_IN_RAM, shellQuote: false}
      label: Max records in ram
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: '500000'
      type: ['null', int]
    - doc: Validation stringency for all sam files read by this program. Setting stringency
        to silent can improve performance when processing a bam file in which variable-length
        data (read, qualities, tags) do not otherwise need to be decoded.
      id: validation_stringency
      inputBinding: {position: 4, prefix: --VALIDATION_STRINGENCY, shellQuote: false}
      label: Validation stringency
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: STRICT
      type:
      - 'null'
      - name: validation_stringency
        symbols: [STRICT, LENIENT, SILENT]
        type: enum
    - doc: Memory which will be allocated for execution.
      id: memory_per_job
      label: Memory Per Job
      sbg:category: Execution
      type: ['null', int]
    - doc: Memory overhead which will be allocated for one job.
      id: memory_overhead_per_job
      label: Memory Overhead Per Job
      sbg:category: Execution
      type: ['null', int]
    - doc: "Sort order of output file.   Required. Possible values: {\n                              queryname
        (Sorts according to the readname. This will place read-pairs and other derived\n
        \                             reads (secondary and supplementary) adjacent
        to each other. Note that the readnames are\n                              compared
        lexicographically, even though they may include numbers. In paired reads,
        Read1\n                              sorts before Read2.)\n                              coordinate
        (Sorts primarily according to the SEQ and POS fields of the record. The\n
        \                             sequence will sorted according to the order
        in the sequence dictionary, taken from from\n                              the
        header of the file. Within each reference sequence, the reads are sorted by
        the\n                              position. Unmapped reads whose mates are
        mapped will be placed near their mates. Unmapped\n                              read-pairs
        are placed after all the mapped reads and their mates.)\n                              duplicate
        (Sorts the reads so that duplicates reads are adjacent. Required that the\n
        \                             mate-cigar (MC) tag is present. The resulting
        will be sorted by library, unclipped 5-prime\n                              position,
        orientation, and mate's unclipped 5-prime position.)\n                              }"
      id: sort_order
      inputBinding: {position: 7, prefix: --SORT_ORDER, shellQuote: false}
      sbg:altPrefix: -SO
      sbg:category: Required  Arguments
      type:
        name: sort_order
        symbols: [queryname, coordinate, duplicate]
        type: enum
    - doc: This input allows a user to set the desired CPU requirement when running
        a tool or adding it to a workflow.
      id: cpu_per_job
      label: CPU per job
      sbg:category: Platform Options
      sbg:toolDefaultValue: '1'
      type: ['null', int]
    - doc: Output file format.
      id: output_file_format
      label: Output file format
      sbg:category: Optional Arguments
      sbg:toolDefaultValue: Same as input
      type:
      - 'null'
      - name: output_file_format
        symbols: [bam, sam]
        type: enum
    label: GATK SortSam
    outputs:
    - doc: Sorted BAM or SAM output file.
      id: out_alignments
      label: Sorted BAM/SAM
      outputBinding: {glob: '*am', outputEval: '$(inheritMetadata(self, inputs.in_alignments))'}
      sbg:fileTypes: BAM, SAM
      secondaryFiles: ["${\n   if (inputs.create_index)\n   {\n       return [self.basename
          + \".bai\", self.nameroot + \".bai\"]\n   }\n   else {\n       return [];
          \n   }\n}"]
      type: ['null', File]
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: "${\n    return inputs.cpu_per_job ?
        inputs.cpu_per_job : 1;\n}", ramMin: "${\n    var memory = 4096;\n    if (inputs.memory_per_job)
        \n    {\n        memory = inputs.memory_per_job;\n    }\n    if (inputs.memory_overhead_per_job)\n
        \   {\n        memory += inputs.memory_overhead_per_job;\n    }\n    return
        memory;\n}"}
    - {class: DockerRequirement, dockerPull: 'images.sbgenomics.com/stefan_stojanovic/gatk:4.1.0.0'}
    - class: InitialWorkDirRequirement
      listing: []
    - class: InlineJavascriptRequirement
      expressionLib: ["var updateMetadata = function(file, key, value) {\n    file['metadata'][key]
          = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata)
          {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n
          \   for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n
          \   }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n
          \   var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2
          = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example
          = o2[i]['metadata'];\n        for (var key in example) {\n            if
          (i == 0)\n                commonMetadata[key] = example[key];\n            else
          {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete
          commonMetadata[key]\n                }\n            }\n        }\n    }\n
          \   if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n
          \   } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i]
          = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar
          toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy
          = function(files, key) {\n    var groupedFiles = [];\n    var tempDict =
          {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n
          \       if (value in tempDict)\n            tempDict[value].push(files[i]);\n
          \       else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict)
          {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar
          orderBy = function(files, key, order) {\n    var compareFunction = function(a,
          b) {\n        if (a['metadata'][key].constructor === Number) {\n            return
          a['metadata'][key] - b['metadata'][key];\n        } else {\n            var
          nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n
          \           if (nameA < nameB) {\n                return -1;\n            }\n
          \           if (nameA > nameB) {\n                return 1;\n            }\n
          \           return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n
          \   if (order == undefined || order == \"asc\")\n        return files;\n
          \   else\n        return files.reverse();\n};", "\nvar setMetadata = function(file,
          metadata) {\n    if (!('metadata' in file))\n        file['metadata'] =
          metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key]
          = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata
          = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2))
          {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n
          \       var example = o2[i]['metadata'];\n        for (var key in example)
          {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n
          \           else {\n                if (!(commonMetadata[key] == example[key]))
          {\n                    delete commonMetadata[key]\n                }\n            }\n
          \       }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1,
          commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++)
          {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n
          \   return o1;\n};"]
    sbg:appVersion: [v1.0]
    sbg:categories: [Utilities, BAM Processing]
    sbg:content_hash: a4d21247730823bddd1b0c24a25cc7b27bea6e061eacc901c23e642f333f458d5
    sbg:contributors: [nens, uros_sipetic]
    sbg:copyOf: veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/19
    sbg:createdBy: uros_sipetic
    sbg:createdOn: 1555498331
    sbg:id: uros_sipetic/gatk-4-1-0-0-demo/gatk-sortsam-4-1-0-0/8
    sbg:image_url: null
    sbg:latestRevision: 8
    sbg:license: Open source BSD (3-clause) license
    sbg:links:
    - {id: 'https://software.broadinstitute.org/gatk/', label: Homepage}
    - {id: 'https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/picard_sam_SortSam.php',
      label: Documentation}
    - {id: 'https://www.ncbi.nlm.nih.gov/pubmed?term=20644199', label: Publications}
    - {id: 'https://github.com/broadinstitute/gatk/', label: Source code}
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1561632457
    sbg:project: uros_sipetic/gatk-4-1-0-0-demo
    sbg:projectName: GATK 4.1.0.0 - Demo
    sbg:publisher: sbg
    sbg:revision: 8
    sbg:revisionNotes: Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/19
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1555498331, 'sbg:revision': 0,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/2}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1555582270, 'sbg:revision': 1,
      'sbg:revisionNotes': Copy of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/9}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557417459, 'sbg:revision': 2, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/11}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557734528, 'sbg:revision': 3, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/13}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558000570, 'sbg:revision': 4, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/14}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558009951, 'sbg:revision': 5, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/15}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558351565, 'sbg:revision': 6, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/17}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1558449641, 'sbg:revision': 7, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/18}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1561632457, 'sbg:revision': 8, 'sbg:revisionNotes': Copy
        of veliborka_josipovic/gatk-4-1-0-0-toolkit-dev/gatk-sortsam-4-1-0-0/19}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Broad Institute
    sbg:toolkit: GATK
    sbg:toolkitVersion: 4.1.0.0
    sbg:validationErrors: []
  sbg:x: 453.0468444824219
  sbg:y: 161.9835205078125
- id: sbg_lines_to_interval_list
  in:
  - {id: input_tsv, source: gatk_createsequencegroupingtsv_4_1_0_0/sequence_grouping}
  label: SBG Lines to Interval List
  out:
  - {id: out_intervals}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    baseCommand: [python, lines_to_intervals.py]
    class: CommandLineTool
    cwlVersion: v1.0
    doc: 'This tools is used for splitting GATK sequence grouping file into subgroups.


      ### Common Use Cases


      Each subgroup file contains intervals defined on single line in grouping file.
      Grouping file is output of GATKs **CreateSequenceGroupingTSV** script which
      is used in best practice workflows sush as **GATK Best Practice Germline Workflow**.'
    id: bix-demo/sbgtools-demo/sbg-lines-to-interval-list/3
    inputs:
    - doc: This file is output of GATKs CreateSequenceGroupingTSV script.
      id: input_tsv
      inputBinding: {position: 1, shellQuote: false}
      label: Input group file
      sbg:category: Required Arguments
      sbg:fileTypes: TSV, TXT
      type: File
    label: SBG Lines to Interval List
    outputs:
    - doc: GATK Intervals files.
      id: out_intervals
      label: Intervals
      sbg:fileTypes: INTERVALS, BED
      type: {items: File, type: array}
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: 1, ramMin: 1000}
    - {class: DockerRequirement, dockerPull: 'python:2.7'}
    - class: InitialWorkDirRequirement
      listing:
      - {entry: "import sys\nimport hashlib\nimport os\nimport json\n\nobj_template
          = {\n    'basename': '',\n    'checksum': '',\n    'class': 'File',\n    'dirname':
          '',\n    'location': '',\n    'nameext': 'intervals',\n    'nameroot': '',\n
          \   'path': '',\n    'size': '',\n}\n\nwith open(sys.argv[1], 'r') as f:\n\n
          \   obj_list = []\n    sys.stderr.write('Reading file {}\\n'.format(sys.argv[1]))\n
          \   nameroot = '.'.join(sys.argv[1].split('/')[-1].split('.')[:-1])\n    for
          i, line in enumerate(f):\n        out_file_name = '{}.group.{}.intervals'.format(nameroot,
          i+1)\n        out_file = open(out_file_name, 'a')\n        for interval
          in line.split():\n            out_file.write(interval + '\\n')\n        out_file.close()\n
          \       sys.stderr.write('Finished writing to file {}\\n'.format(out_file_name))\n\n
          \       obj = dict(obj_template)\n        obj['basename'] = out_file_name\n
          \       obj['checksum'] = 'sha1$' + hashlib.sha1(open(out_file_name, 'r').read()).hexdigest()\n
          \       obj['dirname'] = os.getcwd()\n        obj['location'] = '/'.join([os.getcwd(),
          out_file_name])\n        obj['nameroot'] = '.'.join(out_file_name.split('.')[:-1])\n
          \       obj['path'] = '/'.join([os.getcwd(), out_file_name])\n        obj['size']
          = os.path.getsize('/'.join([os.getcwd(), out_file_name]))\n\n        obj_list.append(obj)\n\n
          \   out_json = {'out_intervals': obj_list}\n\n    json.dump(out_json, open('cwl.output.json',
          'w'), indent=1)\n    sys.stderr.write('Job done.\\n')\n", entryname: lines_to_intervals.py,
        writable: false}
    - {class: InlineJavascriptRequirement}
    sbg:appVersion: [v1.0]
    sbg:content_hash: af381134a8f7ec0eb913f8b421e9ce252619ccbe28180f09bae88646c225d7955
    sbg:contributors: [stefan_stojanovic, nens, uros_sipetic]
    sbg:createdBy: stefan_stojanovic
    sbg:createdOn: 1555578335
    sbg:id: bix-demo/sbgtools-demo/sbg-lines-to-interval-list/3
    sbg:image_url: null
    sbg:latestRevision: 3
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1557936567
    sbg:project: bix-demo/sbgtools-demo
    sbg:projectName: SBGTools - Demo New
    sbg:publisher: sbg
    sbg:revision: 3
    sbg:revisionNotes: output required
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': stefan_stojanovic, 'sbg:modifiedOn': 1555578335, 'sbg:revision': 0,
      'sbg:revisionNotes': null}
    - {'sbg:modifiedBy': stefan_stojanovic, 'sbg:modifiedOn': 1555578408, 'sbg:revision': 1,
      'sbg:revisionNotes': ''}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1556028074, 'sbg:revision': 2,
      'sbg:revisionNotes': Add BED file type}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557936567, 'sbg:revision': 3, 'sbg:revisionNotes': output
        required}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Stefan Stojanovic
    sbg:toolkit: SBG Tools
    sbg:toolkitVersion: '1.0'
    sbg:validationErrors: []
  sbg:x: 995.1798706054688
  sbg:y: 109.42027282714844
- id: sbg_lines_to_interval_list_1
  in:
  - {id: input_tsv, source: gatk_createsequencegroupingtsv_4_1_0_0/sequence_grouping_with_unmapped}
  label: SBG Lines to Interval List
  out:
  - {id: out_intervals}
  run:
    $namespaces: {sbg: 'https://sevenbridges.com'}
    baseCommand: [python, lines_to_intervals.py]
    class: CommandLineTool
    cwlVersion: v1.0
    doc: 'This tools is used for splitting GATK sequence grouping file into subgroups.


      ### Common Use Cases


      Each subgroup file contains intervals defined on single line in grouping file.
      Grouping file is output of GATKs **CreateSequenceGroupingTSV** script which
      is used in best practice workflows sush as **GATK Best Practice Germline Workflow**.'
    id: bix-demo/sbgtools-demo/sbg-lines-to-interval-list/3
    inputs:
    - doc: This file is output of GATKs CreateSequenceGroupingTSV script.
      id: input_tsv
      inputBinding: {position: 1, shellQuote: false}
      label: Input group file
      sbg:category: Required Arguments
      sbg:fileTypes: TSV, TXT
      type: File
    label: SBG Lines to Interval List
    outputs:
    - doc: GATK Intervals files.
      id: out_intervals
      label: Intervals
      sbg:fileTypes: INTERVALS, BED
      type: {items: File, type: array}
    requirements:
    - {class: ShellCommandRequirement}
    - {class: ResourceRequirement, coresMin: 1, ramMin: 1000}
    - {class: DockerRequirement, dockerPull: 'python:2.7'}
    - class: InitialWorkDirRequirement
      listing:
      - {entry: "import sys\nimport hashlib\nimport os\nimport json\n\nobj_template
          = {\n    'basename': '',\n    'checksum': '',\n    'class': 'File',\n    'dirname':
          '',\n    'location': '',\n    'nameext': 'intervals',\n    'nameroot': '',\n
          \   'path': '',\n    'size': '',\n}\n\nwith open(sys.argv[1], 'r') as f:\n\n
          \   obj_list = []\n    sys.stderr.write('Reading file {}\\n'.format(sys.argv[1]))\n
          \   nameroot = '.'.join(sys.argv[1].split('/')[-1].split('.')[:-1])\n    for
          i, line in enumerate(f):\n        out_file_name = '{}.group.{}.intervals'.format(nameroot,
          i+1)\n        out_file = open(out_file_name, 'a')\n        for interval
          in line.split():\n            out_file.write(interval + '\\n')\n        out_file.close()\n
          \       sys.stderr.write('Finished writing to file {}\\n'.format(out_file_name))\n\n
          \       obj = dict(obj_template)\n        obj['basename'] = out_file_name\n
          \       obj['checksum'] = 'sha1$' + hashlib.sha1(open(out_file_name, 'r').read()).hexdigest()\n
          \       obj['dirname'] = os.getcwd()\n        obj['location'] = '/'.join([os.getcwd(),
          out_file_name])\n        obj['nameroot'] = '.'.join(out_file_name.split('.')[:-1])\n
          \       obj['path'] = '/'.join([os.getcwd(), out_file_name])\n        obj['size']
          = os.path.getsize('/'.join([os.getcwd(), out_file_name]))\n\n        obj_list.append(obj)\n\n
          \   out_json = {'out_intervals': obj_list}\n\n    json.dump(out_json, open('cwl.output.json',
          'w'), indent=1)\n    sys.stderr.write('Job done.\\n')\n", entryname: lines_to_intervals.py,
        writable: false}
    - {class: InlineJavascriptRequirement}
    sbg:appVersion: [v1.0]
    sbg:content_hash: af381134a8f7ec0eb913f8b421e9ce252619ccbe28180f09bae88646c225d7955
    sbg:contributors: [stefan_stojanovic, nens, uros_sipetic]
    sbg:createdBy: stefan_stojanovic
    sbg:createdOn: 1555578335
    sbg:id: bix-demo/sbgtools-demo/sbg-lines-to-interval-list/3
    sbg:image_url: null
    sbg:latestRevision: 3
    sbg:modifiedBy: nens
    sbg:modifiedOn: 1557936567
    sbg:project: bix-demo/sbgtools-demo
    sbg:projectName: SBGTools - Demo New
    sbg:publisher: sbg
    sbg:revision: 3
    sbg:revisionNotes: output required
    sbg:revisionsInfo:
    - {'sbg:modifiedBy': stefan_stojanovic, 'sbg:modifiedOn': 1555578335, 'sbg:revision': 0,
      'sbg:revisionNotes': null}
    - {'sbg:modifiedBy': stefan_stojanovic, 'sbg:modifiedOn': 1555578408, 'sbg:revision': 1,
      'sbg:revisionNotes': ''}
    - {'sbg:modifiedBy': uros_sipetic, 'sbg:modifiedOn': 1556028074, 'sbg:revision': 2,
      'sbg:revisionNotes': Add BED file type}
    - {'sbg:modifiedBy': nens, 'sbg:modifiedOn': 1557936567, 'sbg:revision': 3, 'sbg:revisionNotes': output
        required}
    sbg:sbgMaintained: false
    sbg:toolAuthor: Stefan Stojanovic
    sbg:toolkit: SBG Tools
    sbg:toolkitVersion: '1.0'
    sbg:validationErrors: []
  sbg:x: 990.2066040039062
  sbg:y: -57.584449768066406
