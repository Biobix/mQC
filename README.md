<p style="right">
  <img src="https://github.com/Biobix/mQC/blob/master/mqc_tools/logo_mqc2_whitebg.png" width="70" title="mQC logo">
</p>

![mQC logo alt >](https://github.com/Biobix/mQC/blob/master/mqc_tools/logo_mqc2_whitebg.png)

img[alt$=">"] {
  float:right;
}

# MappingQC  (stand-alone version)

MappingQC is a tool to easily generate some figures which give a nice overview of the quality of the mapping of ribosome profiling data. More specific, it gives an overview of the P site offset calculation, the gene distribution and the metagenic classification. Furthermore, MappingQC does a thorough analysis of the triplet periodicity and the linked triplet phase (typical for ribosome profiling) in the canonical transcript of your data. Especially, the link between the phase distribution and the RPF length, the relative sequence position and the triplet identity are taken into account.

## Example

Main script is mQC.pl and you can run it as in following example:
```
mQC.pl --experiment_name yourexperimentname --samfile yoursamfile.sam --cores 20 --species human --ens_v 86 --ens_db ENS_hsa_86.db --unique N --offset plastid --plastid_bam yourbamfile.bam  --tool_dir mqc_tools
```

The stand-alone version of mappingQC (mQC) can be installed from Conda (Be sure to install it in an Python 2.7 environment, not Python3!):
```
conda install mqc
```

More information on Conda can be found on: https://conda.io/docs/user-guide/index.html

A Galaxy version of mQC is also available in the Galaxy tool shed: https://toolshed.g2.bx.psu.edu/view/rnateam/mqc/523dd2a71270. An example mapping procedure in Galaxy, consisting of FastQC on the raw data, next the alignment of the data with STAR and subsequent FastQC and mQC on the aligned data, can be seen in http://galaxy.ugent.be/u/stevenv/h/gao2014hekmqcdemo.

MappingQC is also built in in the PROTEOFORMER pipeline, an automated pipeline for analysis of ribosome profiling sequencing data: https://github.com/Biobix/proteoformer

## Side notes on Conda installation

### Set up channels

Before installing mQC, it is important that the installation channels are added and put in the right order. To do so, it is important to add them in the following order:
```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

### Readline package

The R package of conda-forge can give errors in combination with the readline package version 6.2 of the defaults channel. To solve these errors, just update the readline package:
```
conda update readline
```

## Input parameters

  * work_dir: working directory to run the scripts in (default: current working directory)
  * experiment_name: customly chosen experiment name for the mappingQC run (mandatory)
  * samfile: path to the SAM/BAM file that comes out of the mapping script of PROTEOFORMER (mandatory)
  * cores: the amount of cores to run the script on (integer, default: 5)
  * species: the studied species (mandatory)
    Already implemented species:
      - Human (Homo sapiens)
      - Mouse (Mus musculus)
      - Fruitfly (drosophila melanogaster)
      - Zebrafish (Danio rerio)
      - Yeast (Saccharomyces cerevisiae)
      - Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344
      - c.elegans (Caenorhabditis_elegans)
  * ens_v: the version of the Ensembl database you want to use
  * tmp: temporary folder for storing temporary files of mappingQC (default: work_dir/tmp)
  * unique: whether to use only the unique alignments.
   Possible options: Y, N (default Y)
  * mapper: the mapper you used to generate the SAM file (STAR/TopHat2/HiSat2, default: STAR)
  * maxmultimap: the maximum amount of multimapped positions used for filtering the reads (default: 16)
  * ens_db: path to the Ensembl SQLite database with annotation info. If you want mappingQC to download the right Ensembl database automatically for you, put in 'get' for this parameter (mandatory)
  * offset: the offset determination method.
    Possible options:		
      - plastid: calculate the offsets with Plastid (Dunn et al. 2016)
        The mapping bam file will be needed for Plastid. You can put the path to the BAM file in with the --plastid_bam argument. If you put in '--plastid_bam convert', then mappingQC converts the SAM file argument to a BAM file and uses this one for Plastid. (default: convert)
        Furthermore, you need to give the minimum and maximum RPRF length for Plastid offset generation in the --min_length_plastid and --max_length_plastid arguments respectively. To assess an ideal minimum and maximum value, you can use the length distribution given in the FastQC output (<https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>) (default values: 22 and 34)
      - standard: use the standard offsets from the paper of Ingolia et al. (2012) (default option)
      - from_file: use offsets from an input file
        The offsets input file should be given in the —offset_file argument
  * min_length_gd: minimum RPF length used for gene distributions and metagenic classification (default: 26).
  * max_length_gd: maximum RPF length used for gene distributions and metagenic classification (default: 34).
  * outfolder: the folder to store the output files (default: work_dir/mQC_output)
  * tool_dir: folder with necessary additional mappingQC tools. More information below in the ‘dependencies’ section. (default: search for the default mQC tool directory location in the active conda environment)
  * plotrpftool: the module that will be used for plotting the RPF-phase figure
   Possible options:
      - grouped2D: use Seaborn to plot a grouped 2D bar chart (default)
      - pyplot3D: use mplot3d to plot a 3D bar chart
        This tool can suffer sometimes from Escher effects, as it tries to plot a 3D plot with the 2D software of pyplot and matplotlib.
      - mayavi: use the mayavi package to plot a 3D bar chart
        This tool only works on local systems with graphical cards.
  * outhtml: custom name for the output HTML file (default: work_dir/mQC_experiment_name.html)
  * outzip: custom name for output ZIP file (default: work_dir/mQC_experiment_name.zip)

## Output

MappingQC makes a folder with different plots. The main output file in this folder however is an HTML file that combines all these figures in a nice overview. To transfer your results to another system or location, mappingQC compresses the results folder also to a zip file, so that all results can be unpacked together on the target location.
Following figures are included in the folder and in the overview HTML file:
  * A table with practical information about the mappingQC analysis (analysis time, input parameters, size of the SAM file,..)
  * Plastid offset information (if applicable)
  * Gene distributions: counts summed over genes and shown in 3 different plots (ranked gene counts, cumulative ranked gene counts and gene count density plot).
  * Metagenic plots: counts summed over different annotations. For the non-coding supergroup, a detail plot is also given.
  * Total phase distribution: shows you how the counts are divided over the three reading frames in the canonical translation products of canonical protein-coding transcripts.
  * RPF-phase distribution: shows you how the counts are divided over the three reading frames, but also separated over the different RPF lengths. This is also based on the canonical translation products of canonical protein-coding transcripts.
  * Phase - relative position distribution: shows you how the counts are divided over the three reading frames on a metagenic scale over all canonical reading frames. Relative position 0 means metagenically at the beginning of the sequences; relative position 1 means at the metagenic end of the sequences.
  * Triplet identity plots: shows you the reading frame distribution separated for all possible codons. This can pick up if there are any reading frame distribution biases for specific codons. The resulting amino acid or Start/Stop signal of each codon is given as well.

## Dependencies

As you can see in the command line, mappingQC relies on a tool directory with some additional tools. These include:
* metagenic_piecharts.R				An R tool to plot the metagenic piecharts in R
* quality_plots.R				An R tool to plot the gene distribution quality plots in R
* mQC.py					A python (Python2) script to plot all the other plots and assemble all the output in an HTML overview file.
* simulate_UTR_for_prokaryotes.py       A script to simulate UTR regions in the genes annotation GTF file. Plastid requires untranslated regions in front of canonical start positions and for prokaryotes, these regions need to be simulated.

MappingQC relies also on SQLite and the sqlite3 command line tool for for fetching annotation information out of its Ensembl database. Furthermore, the Plastid tool (Dunn et al. 2016) should be installed if you want to use it for calculating offsets.

MappingQC relies on following Perl modules which have to be installed on your system:
* DBI
* Getopt::Long
* Parallel::ForkManager
* CWD
* Data::Dumper (for debugging purposes)

Furthermore, mappingQC relies on following Python2 modules which have to be installed on your system:
* getopt
* defaultdict (collections)
* sqlite3
* pandas
* numpy
* matplotlib (including pyplot, colors, cm, gridspec, ticker and mplot3d)
* seaborn

!! For the 3D plot (counts as a function of phase and RPF length), you have to make an adaptation in the Python2 libraries of mplot3d. The default axes3d.py script (that will be installed if you download and install mplot3d, the one that your python2 actually uses!) needs to be replaced by the axes3d.py script you can find at https://github.com/Biobix/proteoformer/tree/master/MappingQC/mqc_tools/site-packages. You also need to delete the axes3d.pyc script!
mplot3d is not able to plot 3D barcharts in non-cubic environments and this adapted script will solve this issue.

## More information

For more information about mappingQC: contact Steven.Verbruggen@UGent.be or Gerben.Menschaert@UGent.be

## Copyright

Copyright (C) 2017 Verbruggen Steven & Menschaert Gerben

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
