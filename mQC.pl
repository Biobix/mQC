#!/usr/bin/env perl

$|=1;

#####################################
##	mQC (MappingQC): ribosome profiling mapping quality control tool
##  Author: S. Verbruggen
##  Supervised by: G. Menschaert
##
##	Copyright (C) 2017 S. Verbruggen & G. Menschaert
##
##	This program is free software: you can redistribute it and/or modify
##	it under the terms of the GNU General Public License as published by
##	the Free Software Foundation, either version 3 of the License, or
##	(at your option) any later version.
##
##	This program is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##	GNU General Public License for more details.
##
##	You should have received a copy of the GNU General Public License
##	along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
## 	For more (contact) information visit https://github.com/Biobix/mQC
#####################################

use strict;
use warnings;
use DBI;
use Data::Dumper;
use Getopt::Long;
use v5.10;
use Parallel::ForkManager;
use Cwd;

##############
##Command-line
##############

# nohup perl ./mQC.pl --experiment_name test --samfile untreat.sam --cores 20 --species mouse --ens_db ENS_mmu_86.db --ens_v 86 --offset plastid > nohup_mappingqc.txt &

my($work_dir,$exp_name,$sam,$original_bam,$cores,$species,$version,$tmpfolder,$unique,$mapper,$maxmultimap,$ens_db,$offset_option,$offset_file,$bam,$tool_dir,$plotrpftool,$min_length_plastid,$max_length_plastid,$min_length_gd,$max_length_gd,$outfolder,$outhtml,$outzip,$galaxy,$galaxysam,$galaxytest);
my $help;


GetOptions(
"work_dir:s" => \$work_dir,                 # The working directory                                         Optional argument (default: CWD)
"experiment_name=s" => \$exp_name,          # The experiment name                                           Mandatory argument
"samfile=s"=>\$sam,                         # The samfile/bamfile to do the analysis on                     Mandatory argument
"cores=i"=>\$cores,                         # The amount of cores to use                                    Optional argument (default: 5)
"species=s"=>\$species,                     # The species                                                   Mandatory argument (mouse, human, fruitfly, zebrafish, yeast, SL1344, c.elegans)
"ens_v=i"=>\$version,                       # The Ensembl version                                           Mandatory argument
"tmp:s"=>\$tmpfolder,                       # The tmp folder                                                Optional argument (default: CWD/tmp)
"unique:s"=>\$unique,                       # Consider only unique reads (Y/N)                              Optional argument (default: Y)
"mapper:s"=>\$mapper,                       # The mapper used to generate the SAM file                      Optional argument (default: STAR)
"maxmultimap=i"=>\$maxmultimap,             # The maximum multimapped positions for parsing                 Optional argument (default: 16)
"ens_db=s"=>\$ens_db,                       # The Ensembl db for annotation                                 Mandatory argument
                                                #If this argument is set to 'get', ENS_db will be downloaded first
"offset:s" =>\$offset_option,               # The offset source for parsing alignments                      Optional argument (default: standard)
"offset_file:s" =>\$offset_file,            # The offsets input file                                        Mandatory if offset option equals 'from_file'
"plastid_bam:s" =>\$bam,                    # The corresponding bam file                                    Optional and only used for plastid when offset option equals 'plastid'
                                            # This file will be needed for plastid offset generation. Default bam file will be 'convert' and the bam file will be converted out of the sam file
"min_length_plastid=i" =>\$min_length_plastid,          # Minimum length for plastid offset generation                  Optional argument (default: 22), only used when offset option equals 'plastid'
"max_length_plastid=i" =>\$max_length_plastid,          # Maximum length for plastid offset gene                        Optional argument (default: 34), only used when offset option equals 'plastid'
"min_length_gd=i" =>\$min_length_gd,        #Minimum length for gene distribution and metagenic classification          Optional argument (default: 26)
"max_length_gd=i" =>\$max_length_gd,        #Maximum length for gene distribution and metagenic classification          Optional argument (default: 34)
"tool_dir:s" => \$tool_dir,                 # The directory with all necessary tools                                    Optional argument (default: conda default installation location)
"plotrpftool:s" => \$plotrpftool,           # The module that will be used for plotting the RPF-phase figure
                                                #grouped2D: use Seaborn to plot a grouped 2D bar chart (default)
                                                #pyplot3D: use mplot3d to plot a 3D bar chart (Suffers sometimes from Escher effects)
                                                #mayavi: use the mayavi package to plot a 3D bar chart (only on systems with graphics cards)
"outfolder:s" => \$outfolder,               # The folder for storing output figures                                     Optional argument (default: workdir/mQC_output)
"outhtml:s" => \$outhtml,                   # The output HTML file                                                      Optional argument (default: workdir/mQC_exp_name.html)
"outzip:s" => \$outzip,                     # The output zip file                                                       Optional argument (default: workdir/mQC_exp_name.zip)
"galaxy:s" => \$galaxy,                     # Run through galaxy or not (Y/N)                                           Optional argument (default: N)
"galaxysam:s" => \$galaxysam,               # Parameter needed for galaxy version                                       Optional argument (default: Y)
"galaxytest:s" => \$galaxytest,             # Parameter needed for galaxy version (to run tests)                        Optional argument (default: N)
"help" => \$help                            # Help text option
);

if ($help){
    print_help_text();
    exit;
}

###########################################################################
#Check all input variable and/or get default values and set extra variables
###########################################################################

my $CWD             = getcwd;

# comment on these
if ($work_dir){
    print "Working directory                                        : $work_dir\n";
} else {
    $work_dir = $CWD;
    print "Working directory                                        : $work_dir\n";
}
my $TMP             = ($ENV{'TMP'}) ? $ENV{'TMP'} : ($tmpfolder) ? $tmpfolder : "$CWD/tmp"; # (1) get the TMP environment variable, (2) get the $tmpfolder variable, (3) get current_working_dir/tmp
print "The following tmpfolder is used                          : $TMP\n";
if ($galaxy){
    if ($galaxy ne 'N' && $galaxy ne 'Y'){
        print "ERROR: galaxy option should be Y or N!\n";
        die;
    }
    if ($galaxysam ne 'N' && $galaxysam ne 'Y'){
        print "ERROR: galaxysam option should be Y or N!\n";
        die;
    }
    if ($galaxytest ne 'N' && $galaxytest ne 'Y'){
        print "ERROR: galaxytest option should be Y or N!\n";
        die;
    } elsif (uc($species) ne "HUMAN"&& uc($species) ne "MOUSE" && uc($species) ne "FRUITFLY" && uc($species) ne "ZEBRAFISH" && uc($species) ne "C.ELEGANS") {
        print "ERROR: galaxy test can only run on human, mouse, fruitfly, c.elegans or zebrafish!\n";
        die;
    }
} else {
    $galaxy = 'N';
    $galaxysam = 'Y';
    $galaxytest = 'N';
}


#Check if tmpfolder exists, if not create it...
if (!-e "$TMP") {
    system ("mkdir ". $TMP);
}
if ($exp_name){
    print "The experiment name                                      : $exp_name\n";
} else {
    print_help_text();
    print "\n\n\n";
    die "ERROR: do not forget the experiment name!\n";
}


#Check the extension of the input file
my $ext = "";
if($sam =~ m/\.([^.]+)$/){
    $ext = $1;
    if ($ext eq "sam"){
        if ($sam){
            print "the input sam file                                       : $sam\n";
        } else {
            die "\nDon't forget to pass the bam/sam file!\n\n";
        }
    } elsif ($ext eq "bam"){
        if ($sam){
            print "the input bam file                                       : $sam\n";
            #Convert input bam file to sam format
            system("samtools view -h ".$sam." > ".$TMP."/input.sam");
            $original_bam = $sam;
            $sam = $TMP."/input.sam";
        } else {
            die "\nDon't forget to pass the bam/sam file!\n\n";
        }
    } elsif ($galaxy eq 'Y' && $galaxysam eq 'Y'){
        $ext="sam";
        print "the input sam file                                       : $sam\n";
    } elsif ($galaxy eq 'Y' && $galaxysam eq 'N'){
        $ext="bam";
        print "the input bam file                                       : $sam\n";
        system("samtools view -h ".$sam." > ".$TMP."/input.sam");
        $original_bam = $sam;
        $sam = $TMP."/input.sam";
    } else {
        die "The input file should be in bam/sam format!\n\n";
    }
} else {
    die "Could not match the extension of the input data file $sam !\n\n";
}

if ($species){
    if (uc($species) eq "HUMAN" || uc($species) eq "MOUSE" || uc($species) eq "FRUITFLY" || uc($species) eq "ZEBRAFISH" || uc($species) eq "YEAST" || uc($species) eq "SL1344" || uc($species) eq "C.ELEGANS"){
        print "Species                                                  : $species\n";
    } else {
        die "ERROR: species should be 'human', 'mouse', 'zebrafish', 'yeast', 'SL1344', 'c.elegans' or 'fruifly'!";
    }
} else {
    die "Do not forget the species argument!";
}
if ($unique){
    if ($unique ne "Y" && $unique ne "N"){
        die "ERROR: unique option should be 'Y' or 'N'!";
    }
    print "Consider unique reads:                                   : $unique\n";
} else {
    $unique = "Y";
    print "Consider unique reads:                                   : $unique\n";
}
if ($mapper){
    if ($mapper ne "STAR" && $mapper ne "TopHat2" && $mapper ne "HiSat2"){
        die "ERROR: mapper should be 'STAR' or 'TopHat2' or 'HiSat2'";
    }
    print "Mapper used to generate SAM file:                        : $mapper\n";
} else {
    $mapper = "STAR";
    print "Mapper used to generate SAM file:                        : $mapper\n";
}

#Conversion for species terminolo
my $spec = (uc($species) eq "MOUSE") ? "Mus_musculus" : (uc($species) eq "HUMAN") ? "Homo_sapiens" : (uc($species) eq "SL1344") ? "SL1344" : uc($species) eq "C.ELEGANS" ? "Caenorhabditis_elegans" : (uc($species) eq "ZEBRAFISH") ? "Danio_rerio" : (uc($species) eq "YEAST") ? "Saccharomyces_cerevisiae" : (uc($species) eq "FRUITFLY") ? "Drosophila_melanogaster" : "";
my $spec_short = (uc($species) eq "MOUSE") ? "mmu" : (uc($species) eq "HUMAN") ? "hsa" : (uc($species) eq "ZEBRAFISH") ? "dre" : (uc($species) eq "YEAST") ? "sce" : uc($species) eq "C.ELEGANS" ? "cel" : (uc($species) eq "FRUITFLY") ? "dme" : (uc($species) eq "SL1344") ? "sl1344" : "";
#Old mouse assembly = NCBIM37, new one is GRCm38. Old human assembly = GRCh37, the new one is GRCh38
my $assembly = (uc($species) eq "MOUSE" && $version >= 70 ) ? "GRCm38"
: (uc($species) eq "MOUSE" && $version < 70 ) ? "NCBIM37"
: (uc($species) eq "HUMAN" && $version >= 76) ? "GRCh38"
: (uc($species) eq "HUMAN" && $version < 76) ? "GRCh37"
: (uc($species) eq "ZEBRAFISH") ? "GRCz10"
: (uc($species) eq "SL1344") ? "ASM21085v2"
: (uc($species) eq "YEAST") ? "R64-1-1"
: (uc($species) eq "C.ELEGANS") ? "WBcel235"
: (uc($species) eq "FRUITFLY" && $version < 79) ? "BDGP5"
: (uc($species) eq "FRUITFLY" && $version >= 79) ? "BDGP6" : "";

#UCSC code
my $ucsc;
if ($assembly eq "GRCh38"){
    $ucsc = "hg38";
} elsif ($assembly eq "GRCh37") {
    $ucsc = "hg19";
} elsif ($assembly eq "GRCm38") {
    $ucsc = "mm10";
} elsif ($assembly eq "NCBIM37") {
    $ucsc = "mm9";
} elsif ($assembly eq "BDGP6") {
    $ucsc = "dm6";
} elsif ($assembly eq "GRCz10") {
    $ucsc = "danRer10";
} elsif ($assembly eq "WBcel235") {
    $ucsc = "ce10";
} elsif ($assembly eq "R64-1-1"){
    $ucsc = "sacCer3";
}

#Folder arguments
if ($tool_dir){
    print "The tool directory is set to                             : $tool_dir\n";
} else {
    #Get the conda environment location
    my $env_file = "envs.txt";
    system("conda info -e > ".$env_file);
    my $conda_env = "";
    open(my $fh, "<", $env_file) or die "Could not open conda environments information file!";
    while(my $line = <$fh>){
        chomp($line);
        if($line =~ m/\*\s*(\/.+)$/){
            $conda_env = $1;
        }
    }
    close($fh);
    system("rm -rf ".$env_file);
    if($conda_env eq "" || $conda_env eq "/bin/mqc_tools/"){
        print "Could not find conda environment for default tool directory allocation!\n";
        die;
    }
    $tool_dir = $conda_env."/bin/mqc_tools/";
    print "The conda mqc tool directory is automatically set to     : $tool_dir\n";
}

#Check the scripts you need
if (!-e $tool_dir."/mQC.py"){
    print "Could not find the python mappingQC plotting script mQC.py!\n";
    die;
}
if (!-e $tool_dir."/metagenic_piecharts.R"){
    print "Could not find the R metagenic distribution plotting script metagenic_piecharts.R!\n";
    die;
}
if (!-e $tool_dir."/quality_plots.R"){
    print "Could not find the R quality plotting script quality_plots.R!\n";
    die;
}
if ($outfolder){
    print "The figure output folder is                              : $outfolder\n";
} else {
    $outfolder = $work_dir."/mQC_output/";
    print "The figure output folder is                              : $outfolder\n";
}
if ($plotrpftool){
    if ($plotrpftool eq "grouped2D" || $plotrpftool eq "pyplot3D" || $plotrpftool eq "mayavi"){
        print "RPF phase plotting tool:                                 : $plotrpftool\n";
        if ($plotrpftool eq "pyplot3D"){
            #Check the installation of the pyplot 3D mod
            print "\n\nChecking the installation of pyplot 3D plotting mod\n";
            if (!-e $tool_dir."/install_pyplot3D_mod.py"){
                print "Could not find the python pyplot 3D installation module!\n";
                die;
            }
            system("python ".$tool_dir."/install_pyplot3D_mod.py");
            print "\n\n";
        }
    } else {
        die "The plotrpftool option should be 'grouped2D', 'pyplot3D' or 'mayavi'!\n";
    }
} else {
    $plotrpftool = "grouped2D";
    print "RPF phase plotting tool:                                 : $plotrpftool\n";
}
if ($outhtml){
    print "The output HTML file is                                  : $outhtml\n";
} else {
    $outhtml = $work_dir."/mQC_".$exp_name.".html";
    print "The output HTML file is                                  : $outhtml\n";
}
if ($outzip){
    print "The output zip file is                                   : $outzip\n";
} else {
    $outzip = $work_dir."/mQC_".$exp_name.".zip";
    print "The output zip file is                                   : $outzip\n";
}

#Ensembl options
if ($version){
    print "Ensembl version                                          : $version\n";
    if(uc($species) eq "SL1344"){
        if($version>36){
            print "Error: latest Ensembl Bacteria version is 36!\n";
            die;
        }
    } else {
        if($version>89){
            print "Error: latest Ensembl version is 89!\n";
            die;
        }
    }
} else {
    die "Do not forget the Ensembl version!";
}
if ($ens_db){
    if ($ens_db eq "get"){
        #Download ensembl db
        print "Download Ensembl DB: ".$species." (version ".$version.")\n";
        system("python ".$tool_dir."/ENS_db.py -v ".$version." -s ".$species);
        $ens_db = "ENS_".$spec_short."_".$version.".db";
        #Move ensembl db to tmp folder
        system("mv ".$ens_db." ".$TMP);
        $ens_db = $TMP."/".$ens_db;
    }
    print "the Ensembl DB                                           : $ens_db\n";
} else {
    die "\nDon't forget to pass the Ensembl DB!\n\n";
}
if ($maxmultimap){
    print "Maximun number of loci for reads to be acceptable        : $maxmultimap\n";
} else {
    $maxmultimap = 16;
    print "Maximun number of loci for reads to be acceptable        : $maxmultimap\n";
}
if ($offset_option) {
    if ($offset_option eq "standard" || $offset_option eq "from_file" || $offset_option eq "plastid") {
        print "Offset source                                            : $offset_option\n";
    } else {
        die "Offset argument needs to be \" standard\", \"from_file\" or \"plastid\"!";
    }
} else {
    $offset_option = "standard";
    print "Offset source                                            : $offset_option\n";
}
if ($offset_option eq "from_file"){
    if ($offset_file) {
        print "Offset input file                                        : $offset_file\n";
    } else {
        die "Do not forget the offset input file if offset argument is \"from_file\"!";
    }
} else {
    $offset_file = "";
}
if ($offset_option eq "plastid"){
    if ($bam){
        if ($bam eq "convert"){
            print "Bam file for plastid offset generation will be converted out of sam file\n";
        } else {
            print "Bam file for plastid offset generation                   : $bam\n";
        }
    } else {
        $bam = "convert";
        print "Bam file for plastid offset generation will be converted out of sam file\n";
    }
    if ($min_length_plastid){
        print "Minimum length for plastid offset generation             : $min_length_plastid\n";
    } else {
        $min_length_plastid = 22;
        print "Minimum length for plastid offset generation             : $min_length_plastid\n";
    }
    if ($max_length_plastid){
        print "Maximum length for plastid offset generation             : $max_length_plastid\n";
    } else {
        $max_length_plastid = 34;
        print "Maximum length for plastid offset generation             : $max_length_plastid\n";
    }
}
if ($min_length_gd){
    print "Minimum length for gene distribution and metagene analysis : $min_length_gd\n";
} else {
    $min_length_gd = 26;
    print "Minimum length for gene distribution and metagene analysis : $min_length_gd\n";
}
if ($max_length_gd){
    print "Maximum length for gene distribution and metagene analysis : $max_length_gd\n";
} else {
    $max_length_gd = 34;
    print "Maximum length for gene distribution and metagene analysis : $max_length_gd\n";
}
if ($cores) {
    print "Number of cores to use for analysis                      : $cores\n";
} else {
    $cores = 5;
    print "Number of cores to use for analysis                      : $cores\n";
}

#Download ChromInfo.txt (cfr. get_igenomes.py script PROTEOFORMER)
if (! -e $TMP."/ChromInfo.txt"){
    download_chrominfo($TMP, $ucsc);
} else {
    print "Chromosomal info file already present\n";
}

# Get chromosomes and correct coord_system_id
print "Get chromosomes and coord_system_id...\n";
my $chromosome_sizes; my $coord_system_id; my @ch;

$chromosome_sizes = $TMP."/ChromInfo.txt";
my %chr_sizes = %{get_chr_sizes($chromosome_sizes)};
$coord_system_id = get_coord_system_id($ens_db,$assembly);

#Test on galaxy should run only on Y chromosome
if($galaxytest eq 'Y'){
    my %chr_sizesY;
    $chr_sizesY{'Y'} = $chr_sizes{'Y'};
    %chr_sizes = %chr_sizesY;
}
    
#Download chromosome sequences
if (! -e $TMP."/Chromosomes"){
    system("mkdir ".$TMP."/Chromosomes");
}
if (! -e $TMP."/Chromosomes_BIN"){
    system("mkdir ".$TMP."/Chromosomes_BIN");
}

my $chrom_dir = $TMP."/Chromosomes";
my $BIN_chrom_dir = $TMP."/Chromosomes_BIN";

print "Get/check chromosome fasta files\n";
my $cores_download; #Max 15
if ($cores>15){
    $cores_download = 15;
} else {
    $cores_download = $cores;
}
my $pm_download = new Parallel::ForkManager($cores_download);
print "   Using ".$cores_download." core(s) (max 15)\n   ----------------------\n";

foreach my $chr (keys %chr_sizes){
    
    ### Start parallel process
    $pm_download->start and next;
    
    ### Download chromosome per process
    if (! -e $TMP."/Chromosomes/".$chr.".fa"){
        downloadChromosomeFasta($chr, $spec, $version, $assembly);
    } else {
        print "\t\t*) Chromosome ".$chr." already present\n";
    }
    
    ### Finish
    $pm_download->finish;
    
}

$pm_download->wait_all_children;
print "\n";


########
# MAIN #
########

# Start time
my $start = time;

## Get chromosomes based on seq_region_id ##
# Sqlite Ensembl
my $db_ENS  = $ens_db;
my $dsn_ENS = "DBI:SQLite:dbname=$db_ENS";
my $us_ENS  = "";
my $pw_ENS  = "";
my $chrs = get_chrs($dsn_ENS,$us_ENS,$pw_ENS,\%chr_sizes,$assembly);

#Test on galaxy should run only on Y chromosome
if ($galaxytest eq 'Y') {
    my $chrsY = {};
    $chrsY->{'Y'}->{'seq_region_id'} = $chrs->{'Y'}->{'seq_region_id'};
    $chrs = $chrsY;
}

# Create binary chromosomes if they don't exist
print "\nChecking/Creating binary chrom files ...\n";
create_BIN_chromosomes($BIN_chrom_dir,$cores,$chrs,$work_dir,$TMP);

#Sam file splitting
print "\n";
if (! -e $TMP."/mappingqc"){
    system("mkdir ".$TMP."/mappingqc");
}
my @splitsam = split(/\//, $sam );
my $samFileName = $splitsam[$#splitsam];
@splitsam = split(/\./,$samFileName);
$samFileName = $splitsam[0];
my $samfilechr1 = $TMP."/mappingqc/".$samFileName."_1.sam";

if (-e $samfilechr1){
    print "Splitted sam files already exist\n";
} else {
    print "Splitting genomic mapping per chromosome\n";
    split_SAM_per_chr(\%chr_sizes,$work_dir,$sam, $unique, $mapper);
}

# Construct p offset hash
print "\n";
my $offset_hash = {};
if($offset_option eq "plastid"){
    
    $offset_hash = run_plastid($bam, $TMP, $version, $spec, $assembly, $exp_name, $min_length_plastid, $max_length_plastid);
    
} elsif($offset_option eq "from_file"){
    #Init
    $offset_hash->{"min"} = 1000;
    $offset_hash->{"max"} = 0;
    #Read in file
    open(my $FR, $offset_file) or die "Could not open $offset_file";
    
    #Parse
    while(my $line = <$FR>){
        if($line =~ /^(\d+)\s+(\d+)$/){
            my $length = $1;
            my $offset = $2;
            $offset_hash->{$length} = $offset;
            if($length<$offset_hash->{"min"}){
                $offset_hash->{"min"} = $length;
            }
            if($length>$offset_hash->{"max"}){
                $offset_hash->{"max"} = $length;
            }
        }
    }
} else {
    #Standard P site offset options from Ingolia paper (2012) (cfr. suppl methods in that paper)
    if(uc($species) eq 'FRUITFLY'){
        $offset_hash->{25} = 12;
    }
    $offset_hash->{26} = 12;
    $offset_hash->{27} = 12;
    $offset_hash->{28} = 12;
    $offset_hash->{29} = 12;
    $offset_hash->{30} = 12;
    $offset_hash->{31} = 13;
    $offset_hash->{32} = 13;
    $offset_hash->{33} = 13;
    $offset_hash->{34} = 14;
    
    #Boundaries
    if(uc($species) eq 'FRUITFLY'){
        $offset_hash->{"min"} = 25;
    } else {
        $offset_hash->{"min"} = 26;
    }
    $offset_hash->{"max"} = 34;
}

#Write offsets to csv for output html file
print "Save offsets to csv\n";
offsets_to_csv($offset_hash, $TMP);

print "\n\n";


if ((!-e $TMP."/mappingqc/rpf_phase.csv") || (!-e $TMP."/mappingqc/pos_table_all.csv") || (!-e $TMP."/mappingqc/total_triplet.csv") || (!-e $TMP."/mappingqc/rankedgenes.png") || (!-e $TMP."/mappingqc/cumulative.png") || (!-e $TMP."/mappingqc/density.png") || (!-e $TMP."/mappingqc/annotation_coding.png") || (!-e $TMP."/mappingqc/annotation_noncoding.png")){

    print "RIBOSOMAL PARSING\n";
    system("mkdir ".$TMP."/counts");
    # Init multi core
    my $pm = new Parallel::ForkManager($cores);
    print "   Using ".$cores." core(s)\n   ---------------\n";

    foreach my $chr (keys %chr_sizes){
        
        ### Start parallel process
        $pm->start and next;
        
        ### RIBO parsing
        RIBO_parsing_genomic_per_chr($work_dir,$sam,$chr,$ens_db,$coord_system_id, $offset_hash, $min_length_gd, $max_length_gd);
        
        ### Finish
        print "* Finished chromosome ".$chr."\n";
        $pm->finish;
    }

    # Finish all subprocesses
    $pm->wait_all_children;
    print "\n\n";





    print "PREPARE DATA FOR PLOTTING MODULES\n";

    ## RPF PHASE TABLE ##
    print "\tRPF phase table\n";
    #Count rpf phase table of all chromosomes together
    my $temp_csv_rpf_phase = $TMP."/mappingqc/rpf_phase.csv";
    system("touch ".$temp_csv_rpf_phase);

    my $rpf_phase = {};
    for (my $i=$offset_hash->{'min'};$i<=$offset_hash->{'max'};$i++){
        for (my $j=0;$j<=2;$j++){
            $rpf_phase->{$i}->{$j} = 0;
        }
    }
    foreach my $chr (keys %chr_sizes){
        my $infile = $TMP."/mappingqc/rpf_phase_".$chr.".csv";
        open(IN,"<".$infile) or die $!;
        while(my $line = <IN>){
            my @linesplit = split(',',$line);
            $rpf_phase->{$linesplit[0]}->{0} = $rpf_phase->{$linesplit[0]}->{0} + $linesplit[1];
            $rpf_phase->{$linesplit[0]}->{1} = $rpf_phase->{$linesplit[0]}->{1} + $linesplit[2];
            $rpf_phase->{$linesplit[0]}->{2} = $rpf_phase->{$linesplit[0]}->{2} + $linesplit[3];
        }
        close(IN);
        system("rm -rf ".$infile);
    }

    #Write rpf phase table to temp csv
    open(OUT_PHASE, "+>>".$temp_csv_rpf_phase);
    foreach my $rpf (keys %{$rpf_phase}){
        print OUT_PHASE $rpf.",".$rpf_phase->{$rpf}->{0}.",".$rpf_phase->{$rpf}->{1}.",".$rpf_phase->{$rpf}->{2}."\n";
    }
    close(OUT_PHASE);

    ## PHASE RELATIVE POSITION DISTRIBUTION
    print "\tPhase - relative position distribution\n";
    #Cat phase-position tmp files
    my $temp_csv_all_pos = $TMP."/mappingqc/pos_table_all.csv";
    system("touch ".$temp_csv_all_pos);

    foreach my $chr (keys %chr_sizes){
        my $temp_csv_chr_pos = $TMP."/mappingqc/phase_position_".$chr.".csv";
        system("cat ".$temp_csv_chr_pos." >> ".$temp_csv_all_pos);
        system("rm -rf ".$temp_csv_chr_pos);
    }

    ## TRIPLET IDENTITY PHASE FILE
    print "\tTriplet identity distributions\n";
    #Read in chr tmp files
    my $triplet_phase = {};
    foreach my $chr (keys %chr_sizes){
        my $infile = $TMP."/mappingqc/triplet_phase_".$chr.".csv";
        open(IN, "< ".$infile) or die $!;
        while(my $line = <IN>){
            chomp($line);
            my @linesplit = split(',',$line);
            my $triplet = $linesplit[0];
            my $phase = $linesplit[1];
            my $count = $linesplit[2];
            if(exists $triplet_phase->{$triplet}->{$phase}){
                $triplet_phase->{$triplet}->{$phase} = $triplet_phase->{$triplet}->{$phase} + $count;
            } else {
                $triplet_phase->{$triplet}->{$phase} = $count;
            }
        }
        close(IN);
        system("rm -rf ".$infile);
    }

    #Write total file for triplet identity
    my $temp_total_triplet = $TMP."/mappingqc/total_triplet.csv";
    open(OUT_TOTAL_TRIPLET, "+>> ".$temp_total_triplet);
    foreach my $triplet (keys %{$triplet_phase}){
        foreach my $phase (keys %{$triplet_phase->{$triplet}}){
            print OUT_TOTAL_TRIPLET $triplet.",".$phase.",".$triplet_phase->{$triplet}->{$phase}."\n";
        }
    }
    close(OUT_TOTAL_TRIPLET);

} else {
    print "Ribosomal parsing already done\n"
}

###########
## Gene distributions
###########

if((!-e $TMP."/mappingqc/rankedgenes.png") || (!-e $TMP."/mappingqc/cumulative.png") || (!-e $TMP."/mappingqc/density.png")){
    print "\nGene distribution\n";
    gene_distribution($db_ENS, $us_ENS, $pw_ENS, \%chr_sizes, $cores, $coord_system_id, $tool_dir);
} else {
    print "\nGene distribution already constructed\n";
}

###########
## Metagenic classification
###########

if((!-e $TMP."/mappingqc/annotation_coding.png") || (!-e $TMP."/mappingqc/annotation_noncoding.png")){
    print "\nMetagenic classification\n";
    metagenic_analysis($db_ENS, $us_ENS, $pw_ENS, \%chr_sizes, $cores, $coord_system_id, $tool_dir);
} else {
    print "\nMetagenic classification already done\n";
}

#Run python plotting script
print "\n\n\n\n";
print "Run python plotting script\n";
if ($ext eq "bam"){
    $sam = $original_bam;
}
my $python_command = "python ".$tool_dir."/mQC.py -g ".$galaxy." -a ".$galaxysam." -y ".$galaxytest." -t ".$TMP." -s ".$sam." -n ".$exp_name." -o ".$outfolder." -h ".$outhtml." -z ".$outzip." -p \"".$offset_option."\" -e ".$ens_db." -d ".$species." -v ".$version." -u ".$unique." -x ".$plotrpftool;
if ($offset_option eq "plastid"){
    my $offset_img = $TMP."/plastid/".$exp_name."_p_offsets.png";
    $python_command = $python_command." -i ".$offset_img;
}
print "Python command:\n\t";
print $python_command."\n";
system($python_command);


# End time
print "   DONE! \n";
my $end = time - $start;
printf("Runtime: %02d:%02d:%02d\n\n",int($end/3600), int(($end % 3600)/60), int($end % 60));



############
# THE SUBS #
############

## Gene distribution for each chromosome ##
sub gene_distribution_chr{
    
    #Catch
    my $db_ens = $_[0];
    my $us_ens = $_[1];
    my $pw_ens = $_[2];
    my $chr = $_[3];
    my $coord_system_id = $_[4];
    
    #Open files
    my $out_chr_table = $TMP."/mappingqc/genedistribution_".$chr.".txt";
    open OUT_CHR_GD,"+>>".$out_chr_table or die $!;
    
    #Connect to ensembl db
    my $dsn_ens = "DBI:SQLite:dbname=".$db_ens;
    my $dbh = dbh($dsn_ens, $us_ens, $pw_ens);
    
    #Get seq region id
    my $seq_region = get_seq_region_id($dbh, $chr, $coord_system_id);
    
    #Get all genes with start and stop position
    my $query1 = "SELECT stable_id,seq_region_start,seq_region_end,seq_region_strand FROM gene WHERE seq_region_id = '$seq_region'";
    my $execute1 = $dbh->prepare($query1);
    $execute1->execute();
    
    my %genes;
    while(my @result1 = $execute1->fetchrow_array()){
        #$result1[0]: gene stable_id
        #$result1[1]: gene seq_region_start
        #$result1[2]: gene seq_region_end
        #$result1[3]: gene seq_region_strand
        
        $genes{$chr.":".$result1[3]}{$result1[0]} = [$result1[1],$result1[2]];
    }
    $execute1->finish();
    
    #Disconnect from ensembl
    $dbh->disconnect();
    
    #Make lists of genes (forward and reverse) sorted based on coordinates
    my %for_genes = %{$genes{$chr.":1"}};
    my %rev_genes = %{$genes{$chr.":-1"}};
    my(@genes_for,@genes_rev);
    foreach my $gene_id (sort { $genes{$chr.":1"}{$a}[0] <=> $genes{$chr.":1"}{$b}[0] } keys(%for_genes)){
        push(@genes_for,$gene_id);
    }
    foreach my $gene_id (sort { $genes{$chr.":-1"}{$a}[0] <=> $genes{$chr.":-1"}{$b}[0] } keys(%rev_genes)){
        push(@genes_rev,$gene_id);
    }
    
    ##############
    ## RIBO-SEQ -> READs (~A-site position): determine gene distribution
    ##############
    print "\t\tAnnotating ribo-seq reads of chr ".$chr."\n";
    
    # Get ribo-seq reads, split per strand
    my ($ribo_for, $pos_for) = get_reads($chr,1);
    my ($ribo_rev, $pos_rev) = get_reads($chr,-1);
    
    #Init
    my %gene_count;
    my $intergenic_count;
    
    # Loop over forward ribo-seq reads
    my @window_genes_for = (); # Init window with genes
    foreach my $pos (@{$pos_for}){
        #Push all genes into window where start<window_pos
        my $last_added_index_g = -1;
        foreach my $gene_id (@genes_for){
            if($genes{$chr.":1"}{$gene_id}[0] <= $pos){
                $last_added_index_g++;
                push(@window_genes_for,$gene_id);
            }else{
                #Don't unnecessarily loop over all genes
                last;
            }
        }
        
        #Get rid of gene ids in gene list already in the window
        splice(@genes_for, 0, $last_added_index_g+1);
        
        #Get rid of gene ids in window where end coordinate < window position
        @window_genes_for = grep {$genes{$chr.":1"}{$_}[1] >= $pos} @window_genes_for;
        
        #Annotate read count for all genes in the window
        my $def = 0;
        foreach my $gene_id (@window_genes_for){
            $def++; #defined in genes
            $gene_count{$gene_id} += $ribo_for->{$pos}{'count'};
        }
        
        if($def==0){
            $intergenic_count += $ribo_for->{$pos}{'count'};
        }
    }#Close forward loop
    
    my @window_genes_rev = ();
    #Loop over reverse ribo-seq reads
    foreach my $pos (@{$pos_rev}){
        #Push all genes into window where start<window_pos
        my $last_added_index_g = -1;
        foreach my $gene_id (@genes_rev){
            if($genes{$chr.":-1"}{$gene_id}[0] <= $pos){
                $last_added_index_g++;
                push(@window_genes_rev,$gene_id);
            }else{
                #Don't unnecessarily loop over all genes
                last;
            }
        }
        
        #Get rid of gene ids in gene list already in the window
        splice(@genes_rev, 0, $last_added_index_g+1);
        
        #Get rid of gene ids in window where end coordinate < window position
        @window_genes_rev = grep {$genes{$chr.":-1"}{$_}[1] >= $pos} @window_genes_rev;
        
        #Annotate read count for all genes in the window
        my $def = 0;
        foreach my $gene_id (@window_genes_rev){
            $def++; #defined in genes
            $gene_count{$gene_id} += $ribo_rev->{$pos}{'count'};
        }
        
        if($def==0){
            $intergenic_count += $ribo_rev->{$pos}{'count'};
        }
    }#Close reverse loop
    
    ##############
    ## RESULTS: Make table
    ##############
    
    foreach my $id (keys(%gene_count)){
        print OUT_CHR_GD $id."\t".$gene_count{$id}."\n";
    }
    close(OUT_CHR_GD);
    
    print "\t*) Finished gene distribution construction for chromosome ".$chr."\n";
    
    return;
}

## Gene distribution ##
sub gene_distribution{
    
    #Catch
    my $db_ens = $_[0];
    my $us_ens = $_[1];
    my $pw_ens = $_[2];
    my %chr_sizes = %{$_[3]};
    my $cores = $_[4];
    my $coord_system_id = $_[5];
    my $tool_dir = $_[6];
    
    # Open files
    my $out_table = $TMP."/mappingqc/genedistribution.txt";
    system("rm -rf ".$out_table);
    system("touch ".$out_table);
    open(OUT_GD,"+>>".$out_table);
    print OUT_GD "GeneID\tread_count\n";
    close(OUT_GD);
    
    #Init multi core
    my $processes = $cores; # Nr of processes
    my $pm = new Parallel::ForkManager($processes); # Open fork
    
    # Loop through chromosomes
    foreach my $chr(keys(%chr_sizes)){
        print "\tStarting analysis for chr ".$chr."...\n";
        # Start parallel process
        $pm->start and next;
        
        #Chromosomal gene distribution construction
        gene_distribution_chr($db_ens, $us_ens, $pw_ens, $chr, $coord_system_id);
        
        # Close process
        $pm->finish;
    }
    
    # Finish forking
    $pm->wait_all_children;
    
    #Concatenate all chromosomal out tables
    foreach my $chr(keys(%chr_sizes)){
        my $chr_out_table = $TMP."/mappingqc/genedistribution_".$chr.".txt";
        my $tmp_file = $TMP."/mappingqc/tmp.txt";
        system("cat ".$out_table." ".$chr_out_table." > ".$tmp_file);
        system("mv ".$tmp_file." ".$out_table);
        system("rm -rf ".$chr_out_table);
    }
    
    #Make plots
    print "Make gene distribution plots\n";
    my $out_png1 = $TMP."/mappingqc/rankedgenes.png";
    my $out_png2 = $TMP."/mappingqc/cumulative.png";
    my $out_png3 = $TMP."/mappingqc/density.png";
    make_plots_gd($out_table, $out_png1, $out_png2, $out_png3, $tool_dir);
    
    return;
}

## Make plots in R ##
sub make_plots_gd{
    
    #Catch
    my $out_table = $_[0];
    my $out_png1 = $_[1];
    my $out_png2 = $_[2];
    my $out_png3 = $_[3];
    my $tool_dir = $_[4];
    
    #Execute R script
    system("Rscript ".$tool_dir."/quality_plots.R ".$out_table." ".$out_png1." ".$out_png2." ".$out_png3);
    
    return;
}

## Metagenic analysis: chromosomal ##
sub metagenic_analysis_chr{
    
    #Catch
    my $db_ens = $_[0];
    my $us_ens = $_[1];
    my $pw_ens = $_[2];
    my $chr = $_[3];
    my $coord_system_id = $_[4];
    
    #####
    # Ensembl info
    ####
    
    #Connect to ensembl
    my $dbh = DBI->connect('DBI:SQLite:'.$db_ens,$us_ens,$pw_ens,
                    {RaiseError => 1},) || die "Database connection not made: $DBI::errstr";
    
    #chr nomenclature ens
    if (uc($species) eq "FRUITFLY"){
        if($chr eq "M"){
            $chr = "dmel_mitochondrion_genome";
        }
    } elsif (uc($species) eq "YEAST") {
        if($chr eq "MT"){
            $chr = "Mito";
        }
    } elsif (uc($species) eq "ZEBRAFISH"){
        if($chr eq "MT"){
            $chr = "MtDNA";
        }
    }
    
    #Get seq_region_id
    my $query = "SELECT seq_region_id FROM seq_region WHERE coord_system_id = '$coord_system_id' AND name = '$chr'";
    my $execute = $dbh->prepare($query);
    $execute->execute();
    my $seq_region;
    while(my @result = $execute->fetchrow_array()){
        $seq_region = $result[0];
    }
    $execute->finish();
    
    #change back chr nomenclature ens
    if (uc($species) eq "FRUITFLY"){
        if($chr eq "dmel_mitochondrion_genome"){
            $chr = "M";
        }
    } elsif (uc($species) eq "YEAST") {
        if($chr eq "Mito"){
            $chr = "MT";
        }
    } elsif (uc($species) eq "ZEBRAFISH"){
        if($chr eq "MtDNA"){
            $chr = "MT";
        }
    }
    
    #######
    # Transcripts
    #######
    
    #Init
    my $trs_c = {};
    my $trs_nc = {};
    
    #Get coding transcripts
    my $query1 = "SELECT transcript_id,gene_id,seq_region_start,seq_region_end,seq_region_strand,biotype,stable_id FROM transcript WHERE seq_region_id = '$seq_region' AND biotype = 'protein_coding'";
    my $execute1 = $dbh->prepare($query1);
    $execute1->execute();
    
    $trs_c = $execute1->fetchall_hashref('transcript_id');
    $execute1->finish();
    
    # Get all other transcripts
    my $query2 = "SELECT transcript_id,gene_id,seq_region_start,seq_region_end,seq_region_strand,biotype,stable_id FROM transcript WHERE seq_region_id = '$seq_region' AND biotype NOT LIKE '%protein_coding%'";
    my $execute2 = $dbh->prepare($query2);
    $execute2->execute();
    
    $trs_nc = $execute2->fetchall_hashref('transcript_id');
    $execute2->finish();
    
    # Split transcripts in forward and reverse arrays
    my ($tr_c_for,$tr_c_rev);
    foreach my $tr_id (sort { $trs_c->{$a}{'seq_region_start'} <=> $trs_c->{$b}{'seq_region_start'} } keys %{$trs_c}){
        if($trs_c->{$tr_id}{'seq_region_strand'} eq '1'){	# Forward
            push(@$tr_c_for,$tr_id);
        }else{	# Reverse
            push(@$tr_c_rev,$tr_id);
        }
    }
    my ($tr_nc_for,$tr_nc_rev);
    foreach my $tr_id (sort { $trs_nc->{$a}{'seq_region_start'} <=> $trs_nc->{$b}{'seq_region_start'} } keys %{$trs_nc}){
        if($trs_nc->{$tr_id}{'seq_region_strand'} eq '1'){	# Forward strand
            push(@$tr_nc_for,$tr_id);
        }else{	# Reverse strand
            push(@$tr_nc_rev,$tr_id);
        }
    }
    
    ###########
    ## EXONs & UTRs
    ###########
    ## PROTEIN-CODING TRANSCRIPTS
    ###########
    # Init
    my $exon = {};
    my (%cds_for,%cds_rev); # CDS = UTRs + EXONs
    
    foreach my $tr_id(keys %{$trs_c}){
        # Get all exons for obtained transcripts from tables exon and exon_transcript
        my $query3 = "SELECT a.exon_id,b.seq_region_start,b.seq_region_end,b.seq_region_strand,a.rank FROM exon_transcript a JOIN exon b ON a.exon_id = b.exon_id WHERE a.transcript_id = '$tr_id'";
        my $execute3 = $dbh->prepare($query3);
        $execute3->execute();
        $exon = $execute3->fetchall_hashref('exon_id');
        
        $execute3->execute();
        my $highest_rank=0;
        while(my @result3 = $execute3->fetchrow_array()){
            #$result3[0]: exon_id
            #$result3[1]: seq_region_start
            #$result3[2]: seq_region_end
            #$result3[3]: seq_region_strand
            #$result3[4]: rank
            
            for(my $i1=$result3[1];$i1<=$result3[2];$i1++){
                if($result3[3] == 1){	# Forward strand
                    $cds_for{$chr.':'.$i1}{'exon'}++;
                }else{	# Reverse strand
                    $cds_rev{$chr.':'.$i1}{'exon'}++;
                }
            }
            
            #Search highest rank
            if($result3[4]>$highest_rank){
                $highest_rank = $result3[4];
            }
        }
        $execute3->finish();
        
        # Get first and last exon of each transcript from table translation
        my $query4 = "SELECT start_exon_id,end_exon_id,seq_start,seq_end FROM translation WHERE transcript_id = '$tr_id'";
        my $execute4 = $dbh->prepare($query4);
        $execute4->execute();
        
        while(my @result4 = $execute4->fetchrow_array()){
            #$result4[0]: start_exon_id
            #$result4[1]: end_exon_id
            #$result4[2]: seq_start (offset position translation start)
            #$result4[3]: seq_end (offset position translation stop)
            
            my ($start_id,$end_id,$seq_start,$seq_end) = ($result4[0],$result4[1],$result4[2],$result4[3]);
            my ($start_first_exon,$stop_first_exon,$start_last_exon,$stop_last_exon,$rank_first_exon,$rank_last_exon);
            my ($start_codon,$stop_codon);
            
            # Determine start codon and stop codon
            # Determine 5' and 3' UTRs
            
            if($trs_c->{$tr_id}{'seq_region_strand'} eq '1'){ # Forward strand
                $start_first_exon = $exon->{$start_id}{'seq_region_start'};
                $stop_first_exon = $exon->{$start_id}{'seq_region_end'};
                $rank_first_exon = $exon->{$start_id}{'rank'};
                $start_last_exon = $exon->{$end_id}{'seq_region_start'};
                $stop_last_exon = $exon->{$end_id}{'seq_region_end'};
                $rank_last_exon = $exon->{$end_id}{'rank'};
                $start_codon = $start_first_exon + $seq_start - 1;
                $stop_codon = $start_last_exon + $seq_end - 1;
                
                for(my $l1=$start_first_exon;$l1<$start_codon;$l1++){
                    $cds_for{$chr.':'.$l1}{'5UTR'}++;
                }
                #5UTR also in exons before first translated exon
                if ($rank_first_exon>1){
                    foreach my $exon_id (keys(%{$exon})){
                        if ($exon->{$exon_id}->{'rank'} < $rank_first_exon){
                            for(my $l1a = $exon->{$exon_id}->{'seq_region_start'};$l1a<=$exon->{$exon_id}->{'seq_region_end'};$l1a++){
                                $cds_for{$chr.':'.$l1a}{'5UTR'}++;
                            }
                        }
                    }
                }
                
                #3UTR also in exons after the last translated exon
                for(my $l2=($stop_codon+1);$l2<=$stop_last_exon;$l2++){
                    $cds_for{$chr.':'.$l2}{'3UTR'}++;
                }
                if ($rank_last_exon < $highest_rank){
                    foreach my $exon_id (keys(%{$exon})){
                        if ($exon->{$exon_id}->{'rank'} > $rank_last_exon){
                            for (my $l2a = $exon->{$exon_id}->{'seq_region_start'};$l2a<=$exon->{$exon_id}->{'seq_region_end'};$l2a++){
                                $cds_for{$chr.':'.$l2a}{'3UTR'}++;
                            }
                        }
                    }
                }
                
            }elsif($trs_c->{$tr_id}{'seq_region_strand'} eq '-1'){ # Reverse strand
                $start_first_exon = $exon->{$start_id}{'seq_region_end'};
                $stop_first_exon = $exon->{$start_id}{'seq_region_start'};
                $rank_first_exon = $exon->{$start_id}{'rank'};
                $start_last_exon = $exon->{$end_id}{'seq_region_end'};
                $stop_last_exon = $exon->{$end_id}{'seq_region_start'};
                $rank_last_exon = $exon->{$end_id}{'rank'};
                $start_codon = $start_first_exon - $seq_start + 1;
                $stop_codon = $start_last_exon - $seq_end + 1;
                
                #5UTR also in exons before first translated exon
                for(my $l3=($start_codon+1);$l3<=$start_first_exon;$l3++){
                    $cds_rev{$chr.':'.$l3}{'5UTR'}++;
                }
                if ($rank_first_exon>1){
                    foreach my $exon_id (keys(%{$exon})){
                        if ($exon->{$exon_id}->{'rank'} < $rank_first_exon){
                            for (my $l3a=$exon->{$exon_id}->{'seq_region_start'}; $l3a<=$exon->{$exon_id}->{'seq_region_end'};$l3a++){
                                $cds_rev{$chr.':'.$l3a}{'5UTR'}++;
                            }
                        }
                    }
                }
                
                #3UTR also in exons after the last translated exon
                for(my $l4=$stop_last_exon;$l4<$stop_codon;$l4++){
                    $cds_rev{$chr.':'.$l4}{'3UTR'}++;
                }
                
                
                if ($rank_last_exon < $highest_rank){
                    foreach my $exon_id (keys(%{$exon})){
                        if ($exon->{$exon_id}->{'rank'} > $rank_last_exon){
                            for (my $l4a=$exon->{$exon_id}->{'seq_region_start'};$l4a<=$exon->{$exon_id}->{'seq_region_end'};$l4a++){
                                $cds_rev{$chr.':'.$l4a}{'3UTR'}++;
                            }
                        }
                    }
                }
            }
        }
        $execute4->finish();
        
    } #close transcript id loop
    
    ###########
    ## NON-CODING TRANSCRIPTS
    ###########
    # Get biotypes
    my $query5 = "SELECT biotype FROM transcript WHERE biotype NOT LIKE '%protein_coding%' GROUP BY biotype";
    my $execute5 = $dbh->prepare($query5);
    $execute5->execute();
    
    my %biotypes_nc;
    while(my @results5 = $execute5->fetchrow_array()){
        $biotypes_nc{$results5[0]} = 0;
    }
    $execute5->finish();
    # Disconnect
    $dbh->disconnect();
    
    ###########
    ## RIBO-DATA
    ###########
    ###########
    ## ANNOTATE RIBO-SEQ READS
    ###########
    print "\t\tAnnotating ribo-seq reads of chr ".$chr."\n";
    
    # Get ribo-seq reads, split per strand
    my ($ribo_for, $pos_for) = get_reads($chr,1);
    my ($ribo_rev, $pos_rev) = get_reads($chr,-1);
    
    # Init values
    my ($ribo_reads,$ribo_readsnc) = (0,0);
    my ($ribo_exon,$ribo_intron,$ribo_5utr,$ribo_3utr,$ribo_intergenic) = (0,0,0,0,0);
    
    # Loop over forward ribo-seq reads
    my @window_c = (); # Init window with protein-coding transcripts
    my @window_nc = (); # Init window with non protein-coding transcripts
    foreach my $pos (@{$pos_for}){
        #####
        ## PROTEIN-CODING WINDOW
        ######
        # Push all tr_ids to @window_c where tr_start < window_pos
        my $last_added_index_c = -1;
        foreach my $tr_for_id (@$tr_c_for){
            if($trs_c->{$tr_for_id}{'seq_region_start'} <= $pos){
                $last_added_index_c++;
                push(@window_c,$tr_for_id);
            }else{last;} # Don't unnecessarily loop over all transcripts
        }
        
        # Get rid of tr_c_for elements brought into @window_c
        splice(@$tr_c_for, 0, $last_added_index_c+1); #length parameter = last index (- first index (=0) ) + 1
        
        # Get rid of tr_ids in @window_c where tr_end < window_pos
        @window_c = grep {$trs_c->{$_}{'seq_region_end'} >= $pos} @window_c;
        
        #####
        ## NON-CODING WINDOW
        #####
        # Push all tr_ids to @window_nc where tr_start < window_pos
        my $last_added_index_nc = -1;
        foreach my $tr_for_id (@$tr_nc_for){
            if($trs_nc->{$tr_for_id}{'seq_region_start'} <= $pos){
                $last_added_index_nc++;
                push(@window_nc,$tr_for_id);
            }else{last;} # Don't unnecessarily loop over all transcripts
        }
        
        # Get rid of tr_nc_for elements already in @window_nc
        splice(@$tr_nc_for, 0, $last_added_index_nc+1);
        
        # Get rid of tr_ids in @window_nc where tr_end < window_pos
        @window_nc = grep {$trs_nc->{$_}{'seq_region_end'} >= $pos} @window_nc;
        
        #####
        ## ANNOTATE
        #####
        # Loop over windows and check functional annotation of each ribo-seq read
        $ribo_reads = $ribo_reads + $ribo_for->{$pos}{'count'};
        if(@window_c){
            # Annotate reads in PROTEIN-CODING transcripts
            my $count = $ribo_for->{$pos}{'count'};
            
            # Check 5'UTR
            if(defined($cds_for{$chr.':'.$pos}{'5UTR'})){
                $ribo_5utr = $ribo_5utr + $count;
            }
            
            # Check 3'UTR
            elsif(defined($cds_for{$chr.':'.$pos}{'3UTR'})){
                $ribo_3utr = $ribo_3utr + $count;
            }
            
            # Check Exons
            elsif(defined($cds_for{$chr.':'.$pos}{'exon'})){
                $ribo_exon = $ribo_exon + $count;
            }
            
            # If still not defined -> intronic region
            else{
                $ribo_intron = $ribo_intron + $count;
            }
        }elsif(@window_nc){
            # Annotate reads in NON PROTEIN-CODING transcripts
            $ribo_readsnc = $ribo_readsnc + $ribo_for->{$pos}{'count'};
            
            # Define biotype (if #biotypes>0, take a random/first one)
            my @biotypes; my $i = 0;
            foreach my $tr_for_id (@window_nc){
                $biotypes[$i] = $trs_nc->{$tr_for_id}{'biotype'};
                $i++;
            }
            my $random = int(rand(@biotypes));
            my $biotype = $biotypes[$random];
            $biotypes_nc{$biotype} = $biotypes_nc{$biotype} + $ribo_for->{$pos}{'count'};
        }else{
            $ribo_intergenic = $ribo_intergenic + $ribo_for->{$pos}{'count'};
        }
    } #close forward loop
        
    # Loop over reverse ribo-seq reads
    @window_c = (); # Empty window_c
    @window_nc = (); # Empty window_nc
    foreach my $pos (@{$pos_rev}){
        ######
        ## PROTEIN-CODING WINDOW
        ######
        # Push all tr_ids to @window_c where tr_start < window_pos
        my $last_added_index_c = -1;
        foreach my $tr_rev_id (@$tr_c_rev){
            if($trs_c->{$tr_rev_id}{'seq_region_start'} <= $pos){
                $last_added_index_c++;
                push(@window_c,$tr_rev_id);
            }else{last;} # Don't unnecessarily loop over all transcripts
        }
        
        # Get rid of tr_c_rev elements already in @window_c
        splice(@$tr_c_rev, 0, $last_added_index_c+1);
        
        # Get rid of tr_ids in @window_c where tr_end < window_pos
        @window_c = grep {$trs_c->{$_}{'seq_region_end'} >= $pos} @window_c;
        
        ######
        ## NON-CODING WINDOW
        ######
        # Push all tr_ids to @window_nc where tr_start < window_pos
        my $last_added_index_nc = -1;
        foreach my $tr_rev_id (@$tr_nc_rev){
            if($trs_nc->{$tr_rev_id}{'seq_region_start'} <= $pos){
                $last_added_index_nc++;
                push(@window_nc,$tr_rev_id);
            }else{last;} # Don't unnecessarily loop over all transcripts
        }
        
        # Get rid of tr_nc_for elements already in @window_nc
        splice(@$tr_nc_rev, 0, $last_added_index_nc+1);
        
        # Get rid of tr_ids in @window_nc where tr_end < window_pos
        @window_nc = grep {$trs_nc->{$_}{'seq_region_end'} >= $pos} @window_nc;
        
        ######
        ## ANNOTATE
        ######
        # Loop over windows and check functional annotation of each ribo-seq read
        $ribo_reads = $ribo_reads + $ribo_rev->{$pos}{'count'};
        if(@window_c){
            # Annotate reads in PROTEIN-CODING transcripts
            my $count = $ribo_rev->{$pos}{'count'};
            
            # Check 5'UTR
            if(defined($cds_rev{$chr.':'.$pos}{'5UTR'})){
                $ribo_5utr = $ribo_5utr + $count;
            }
            
            # Check 3'UTR
            elsif(defined($cds_rev{$chr.':'.$pos}{'3UTR'})){
                $ribo_3utr = $ribo_3utr + $count;
            }
            
            # Check Exons
            elsif(defined($cds_rev{$chr.':'.$pos}{'exon'})){
                $ribo_exon = $ribo_exon + $count;
            }
            
            # If still not defined -> intronic region
            else{
                $ribo_intron = $ribo_intron + $count;
            }
        }elsif(@window_nc){
            # Annotate reads in NON PROTEIN-CODING transcripts
            $ribo_readsnc = $ribo_readsnc + $ribo_rev->{$pos}{'count'};
            
            # Define biotype (if #biotypes>0, take a random/first one)
            my @biotypes; my $i = 0;
            foreach my $tr_rev_id (@window_nc){
                $biotypes[$i] = $trs_nc->{$tr_rev_id}{'biotype'};
                $i++;
            }
            my $random = int(rand(@biotypes));
            my $biotype = $biotypes[$random];
            $biotypes_nc{$biotype} = $biotypes_nc{$biotype} + $ribo_rev->{$pos}{'count'};
            
        }else{
            $ribo_intergenic = $ribo_intergenic + $ribo_rev->{$pos}{'count'};
        }
    } # Close rev-loop
    
    # Print results
    print OUT1 $chr."\t".$ribo_reads."\t".$ribo_exon."\t".$ribo_5utr."\t".$ribo_3utr."\t".$ribo_intron."\t".$ribo_readsnc."\t".$ribo_intergenic."\n";
    print OUT2 $chr."\t".$ribo_readsnc;
    foreach my $biotype(sort keys %biotypes_nc){
        print OUT2 "\t".$biotypes_nc{$biotype};
    }
    print OUT2 "\n";
    
    print "\t*) Finished metagenic analysis for chromosome ".$chr."\n";
    
    return;
    
}

## Metagenic analysis: total ##
sub metagenic_analysis {
    
    #Catch
    my $db_ens = $_[0];
    my $us_ens = $_[1];
    my $pw_ens = $_[2];
    my %chr_sizes = %{$_[3]};
    my $cores = $_[4];
    my $coord_system_id = $_[5];
    my $tool_dir = $_[6];
    
    #Get Non coding biotypes out of ensembl
    my %biotypes = %{get_nPCbiotypes($db_ens, $us_ens, $pw_ens)};
    
    #Open files
    my $out_table1 = $TMP."/mappingqc/annotation_coding.txt";
    my $out_table2 = $TMP."/mappingqc/annotation_noncoding.txt";
    if(-e $out_table1){
        system("rm -rf ".$out_table1);
    }
    if(-e $out_table2){
        system("rm -rf ".$out_table2);
    }
    open(OUT1,"+>>",$out_table1) or die $!;
    open(OUT2,"+>>",$out_table2) or die $!;
    
    print OUT1 "chr\tribo\texon\t5utr\t3utr\tintron\tnon_protein_coding\tintergenic\n";
    print OUT2 "chr\tnon_protein_coding";
    foreach my $biotype(sort keys %biotypes){
        print OUT2 "\t".$biotype;
    }
    print OUT2 "\n";
    
    #Init multi core
    my $pm = new Parallel::ForkManager($cores);
    
    foreach my $chr (keys(%chr_sizes)){
        
        #Start parallel process for each chromosome
        $pm->start and next;
        
        #main loop
        metagenic_analysis_chr($db_ens, $us_ens, $pw_ens, $chr, $coord_system_id);
        
        #Finish process
        $pm->finish;
        
    }
    
    #Finish all processes
    $pm->wait_all_children;
    
    #Close output stream
    close(OUT1);
    close(OUT2);
    
    #output figures
    my $out_png1 = $TMP."/mappingqc/annotation_coding.png";
    my $out_png2 = $TMP."/mappingqc/annotation_noncoding.png";
    if(-e $out_png1){
        system("rm -rf ".$out_png1);
    }
    if(-e $out_png2){
        system("rm -rf ".$out_png2);
    }
    
    # Make Pie Charts
    piecharts($out_table1,$out_table2,$out_png1,$out_png2,$tool_dir);
    
    return;
}

#Startup R for plotting pie charts
sub piecharts{
    print "Make Pie Charts...\n";
    
    # Catch
    my $out_table1 = $_[0];
    my $out_table2 = $_[1];
    my $out_png1 = $_[2];
    my $out_png2 = $_[3];
    my $tooldir = $_[4];
    
    # Execute Rscript
    system("Rscript ".$tooldir."/metagenic_piecharts.R ".$out_table1." ".$out_table2." ".$out_png1." ".$out_png2);
}

#Get reads out of tmp files
sub get_reads{
    
    # Catch
    my $chr = $_[0];
    my $strand = $_[1];
    
    # Init
    my $ribo_reads = {};
    my $read_keys = [];
    
    #Open chromosomal read tmp file
    my $reads_file = "";
    if($strand eq 1){
        $reads_file = $TMP."/counts/reads_".$chr."_FOR.csv";
    } else {
        $reads_file = $TMP."/counts/reads_".$chr."_REV.csv";
    }
    open(READS,"<",$reads_file);
    while(my $line = <READS>){
        chomp $line;
        my @elements = split(/,/, $line);
        $ribo_reads->{$elements[0]}->{'count'} = $elements[1];
        push(@{$read_keys}, $elements[0]);
    }
    close(READS);
    
    # Return
    return($ribo_reads, $read_keys);
}

### RIBO PARSE PER CHR ###
sub RIBO_parsing_genomic_per_chr {
    
    #Catch
    my $work_dir = $_[0];
    my $sam = $_[1];
    my $chr = $_[2];
    my $ens_db = $_[3];
    my $coord_system_id = $_[4];
    my $offset_hash = $_[5];
    my $min_l_parsing = $_[6];
    my $max_l_parsing = $_[7];
    
    my @splitsam = split(/\//, $sam );
    my $samFileName = $splitsam[$#splitsam];
    @splitsam = split(/\./,$samFileName);
    $samFileName = $splitsam[0];
    
    #Construct phase library
    my ($phase_lib, $triplet_lib) = construct_phase_lib($chr, $ens_db, $coord_system_id);
    
    #Initialize
    my ($genmatchL,$offset,$start,$intron_total,$extra_for_min_strand);
    my $phase_count_RPF = {};
    for (my $i=$offset_hash->{'min'};$i<=$offset_hash->{'max'};$i++){
        for (my $j=0;$j<=2;$j++){
            $phase_count_RPF->{$i}->{$j} = 0;
        }
    }
    my $phase_count_file = $TMP."/mappingqc/rpf_phase_".$chr.".csv";
    my $phase_count_triplet = {};
    my $triplet_count_file = $TMP."/mappingqc/triplet_phase_".$chr.".csv";
    my $chr_sam_file = $TMP."/mappingqc/".$samFileName."_".$chr.".sam";
    my $pos_file = $TMP."/mappingqc/phase_position_".$chr.".csv";
    my $hits_genomic = {};
    
    open (I,"<".$chr_sam_file) || die "Cannot open ".$chr_sam_file."\n";
    open (OUT_POS, "+>>".$pos_file) or die $!;

    
    while(my $line=<I>){
        #Process alignment line
        my @mapping_store = split(/\t/,$line);
        
        #Get strand specifics
        # Sam flag is bitwise. (0x10 SEQ being reverse complemented)
        # 0x10 = 16 in decimal. -> negative strand.
        my $strand = ($mapping_store[1] & 16) ? "-": "+";
        my $CIGAR = $mapping_store[5];
        #Rewrite strand info
        my $strandAlt;
        if($strand eq "+"){
            $strandAlt = 1;
        } else {
            $strandAlt = -1;
        }
        
        #Parse CIGAR to obtain offset,genomic matching length and total covered intronic region before reaching the offset
        ($offset,$genmatchL,$intron_total,$extra_for_min_strand) = parse_RIBO_CIGAR($CIGAR,$strand,$offset_hash);
        #Determine genomic position based on CIGAR string output and mapping position and direction
        $start = ($strand eq "+") ? $mapping_store[3] + $offset + $intron_total : ($strand eq "-") ? $mapping_store[3] - $offset - $intron_total + $extra_for_min_strand -1 : "";
        
        if($genmatchL>=$offset_hash->{"min"} && $genmatchL<=$offset_hash->{"max"}){
            if(exists $phase_lib->{$strandAlt}->{$start}->{"phase"}){
                #Add for RPF-splitted phase distribution
                $phase_count_RPF->{$genmatchL}->{$phase_lib->{$strandAlt}->{$start}->{"phase"}}++;
                #Print to tmp chr phase-position distribbution
                if(exists $phase_lib->{$strandAlt}->{$start}->{"transcriptomic_pos"}){
                    my $read_phase = $phase_lib->{$strandAlt}->{$start}->{"phase"};
                    print OUT_POS $read_phase.",".$phase_lib->{$strandAlt}->{$start}->{"transcriptomic_pos"}."\n";
                }
            }
            if(exists $triplet_lib->{$strandAlt}->{$start}){
                my $read_phase = $phase_lib->{$strandAlt}->{$start}->{"phase"};
                my $read_triplet = $triplet_lib->{$strandAlt}->{$start};
                if(length($read_triplet)==3){
                    if(exists $phase_count_triplet->{$read_triplet}){
                        $phase_count_triplet->{$read_triplet}->{$read_phase}++;
                    } else {
                        #Initialize if triplet is first time seen
                        $phase_count_triplet->{$read_triplet}->{0} = 0;
                        $phase_count_triplet->{$read_triplet}->{1} = 0;
                        $phase_count_triplet->{$read_triplet}->{2} = 0;
                        $phase_count_triplet->{$read_triplet}->{$read_phase}++;
                    }
                }
            }
            
            #Save counts for gene distribution and metagenic classification
            if($genmatchL >= $min_l_parsing && $genmatchL <= $max_l_parsing){
                if(exists $hits_genomic->{$start}->{$strandAlt}){
                    $hits_genomic->{$start}->{$strandAlt}++;
                } else {
                    $hits_genomic->{$start}->{$strandAlt} = 0;
                    $hits_genomic->{$start}->{$strandAlt}++;
                }
            }
        }

    }
    
    #Stop reading out of input files
    close(I);
    
    #Write read counts to files for later use in gene distributions and metagenic classification
    my $chr_for_reads_file = $TMP."/counts/reads_".$chr."_FOR.csv";
    my $chr_rev_reads_file = $TMP."/counts/reads_".$chr."_REV.csv";
    system("rm -rf ".$chr_for_reads_file);
    system("rm -rf ".$chr_rev_reads_file);
    open (CHR_FOR_READS,"+>>".$chr_for_reads_file) or die $!;
    open (CHR_REV_READS,"+>>".$chr_rev_reads_file) or die $!;
    
    foreach my $pos (sort {$a <=> $b} keys %{$hits_genomic}){
        if(exists $hits_genomic->{$pos}->{1}){
            print CHR_FOR_READS $pos.",".$hits_genomic->{$pos}->{1}."\n";
        }
        if(exists $hits_genomic->{$pos}->{-1}){
            print CHR_REV_READS $pos.",".$hits_genomic->{$pos}->{-1}."\n";
        }
    }
    
    #Stop writing reads out
    close(CHR_FOR_READS);
    close(CHR_REV_READS);
    
    #Write to chromosomal rpf phase tmp file
    open (OUT_RPF_PHASE, "+>>".$phase_count_file) or die $!;
    
    for (my $i=$offset_hash->{'min'};$i<=$offset_hash->{'max'};$i++){
        print OUT_RPF_PHASE $i.",".$phase_count_RPF->{$i}->{0}.",".$phase_count_RPF->{$i}->{1}.",".$phase_count_RPF->{$i}->{2}."\n";
    }
    close(OUT_RPF_PHASE);
    
    #Write to chromosomal triplet phase tmp file
    open (OUT_TRIPLET_PHASE, "+>>".$triplet_count_file) or die $!;
    
    foreach my $triplet (keys %{$phase_count_triplet}){
        foreach my $phase (keys %{$phase_count_triplet->{$triplet}}){
            print OUT_TRIPLET_PHASE $triplet.",".$phase.",".$phase_count_triplet->{$triplet}->{$phase}."\n";
        }
    }
    
    close(OUT_TRIPLET_PHASE);
    
    return;
}

#Get non coding biotypes
sub get_nPCbiotypes{
    
    # Catch
    my $db_ensembl = $_[0];
    my $user = $_[1];
    my $pw = $_[2];
    
    # Connect to Ensembl db
    my $dbh = DBI->connect('DBI:SQLite:'.$db_ensembl,$user,$pw,
    {RaiseError => 1},) || die "Database connection not made: $DBI::errstr";
    
    # Query db
    my $query = "SELECT biotype FROM transcript WHERE biotype NOT LIKE '%protein_coding%' GROUP BY biotype";
    my $execute = $dbh->prepare($query);
    $execute->execute();
    
    my %biotypes;
    while(my @results = $execute->fetchrow_array()){
        $biotypes{$results[0]} = 0;
    }
    
    $execute->finish();
    
    # Disconnect
    $dbh->disconnect();
    
    # Return
    return(\%biotypes)
}


#Parse RIBO_CIGARS to obtain offset,genomic read mapping length and total intronic length before offset is reached
sub parse_RIBO_CIGAR {
    
    #Catch
    my $CIGAR = $_[0];
    my $strand = $_[1];
    
    my $CIGAR_SPLIT = splitCigar($CIGAR);
    my $CIGAR_SPLIT_STR = [];
    @$CIGAR_SPLIT_STR = ($strand eq "-") ? reverse @$CIGAR_SPLIT : @$CIGAR_SPLIT;
    my $op_total = @$CIGAR_SPLIT_STR;
    
    my $genmatchL = 0;
    my $op_count = 0;
    #To keep track of total length of genomic + intron (negative strand, reverse position)
    my $extra_for_min_strand = 0;
    #Loop over operation to get total mapping length to calculate the A-site offset
    # and to get total extra length for min_strand (i.e. S(not 5adapt nor 1stTRIM), N (splicing), M, D,I)
    foreach my $operation (@$CIGAR_SPLIT_STR) {
        my $op_length = $operation->[0];
        my $op_type = $operation->[1];
        $op_count++;
        
        if($op_type =~ /^S$/) {
            #Trim leading substitution if only 1 substitution @ 5'
            if ($op_count == 1 && $op_length == 1) {
                next;
            }
            #Clip trailing adaptor substitution
            elsif ($op_count == $op_total) {
                next;
            }
            #Other substitutions are added to RIBO-read genomic-match length
            #And also added to the total matching count
            else {
                $genmatchL = $genmatchL + $op_length;
                $extra_for_min_strand = $extra_for_min_strand + $op_length;
            }
        }
        #Sum matching operations until the offset is reached, then change status to "Y"
        elsif($op_type =~ /^M$/) {
            $genmatchL = $genmatchL + $op_length;
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
        #Insertions elongate the readL and the insertion size is added to the total matching count
        elsif($op_type =~ /^I$/) {
            $genmatchL = $genmatchL + $op_length;
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
        #Splice intronic regions are added to the extra_for_min_strand
        elsif($op_type =~ /^N$/) {
            $extra_for_min_strand = $extra_for_min_strand + $op_length;
        }
    }
    #print "total mapped sequence = $genmatchL\n";
    my $offset = get_offset($genmatchL, $offset_hash);
    my $match_count_total = 0;
    $op_count = 0;
    my $offset_covered = "N";
    my $intron_total = 0;
    #Loop over operations to caculate the total intron length
    foreach my $operation (@$CIGAR_SPLIT_STR) {
        my $op_length = $operation->[0];
        my $op_type = $operation->[1];
        $op_count++;
        #print "$op_type,$op_length\n";
        if($op_type =~ /^S$/) {
            #Trim leading substitution if only 1 substitution @ 5'
            if ($op_count == 1 && $op_length == 1) {
                next;
            }
            #Clip trailing adaptor substitution
            elsif ($op_count == $op_total) {
                next;
            }
            #Other substitutions are added to RIBO-read genomic-match length
            #And also added to the total matching count
            else {
                $match_count_total = $match_count_total + $op_length;
                if ($match_count_total >= $offset) {
                    $offset_covered = "Y";
                    last;
                }
            }
        }
        #Sum matching operations until the offset is reached, then change status to "Y"
        elsif($op_type =~ /^M$/) {
            $match_count_total = $match_count_total + $op_length;
            if ($match_count_total >= $offset) {
                $offset_covered = "Y";
                last;
            }
        }
        #Sum intronic region lengths untill the offset has been covered by the matching operations
        elsif($op_type =~ /^N$/ && $offset_covered eq "N") {
            $intron_total = $intron_total + $op_length;
            
        }
        #Deletion are not counted for the readL
        elsif($op_type =~ /^D$/) {
            next;
            #$genmatchL = $genmatchL - $op_length;
        }
        #Insertions elongate the readL and the insertion size is added to the total matching count
        elsif($op_type =~ /^I$/) {
            $match_count_total = $match_count_total + $op_length;
            if ($match_count_total >= $offset) {
                $offset_covered = "Y";
                last;
            }
        }
        #print "$match_count_total,$offset_covered,$offset\n";
    }
    
    return($offset,$genmatchL,$intron_total,$extra_for_min_strand)
}

#Construct phase lib out of ensembl info
sub construct_phase_lib{
    
    #Catch
    my $chr = $_[0];
    my $eDB = $_[1];
    my $coord_system_id = $_[2];
    
    #Init
    my $phase_lib = {};
    my $triplet_lib = {};
    my $dsn_sqlite_ens = "DBI:SQLite:dbname=$eDB";
    my $us_sqlite_ens  = "";
    my $pw_sqlite_ens  = "";
    
    #Init dbh
    my $dbh_ens = dbh($dsn_sqlite_ens,$us_sqlite_ens,$pw_sqlite_ens);
    
    #Get seq_region_id
    my $seq_region_id = get_seq_region_id($dbh_ens, $chr, $coord_system_id);
    
    #Get transcripts (canonical protein-coding)
    my $transcripts = get_can_transcripts($dbh_ens, $seq_region_id);
    
    #Get exon structure of each transcript
    foreach my $transcript (@$transcripts){
        my($exon_struct,$strand,$max_tr_rank) = get_exon_struct_transcript($dbh_ens, $transcript, $chr);
        ($phase_lib, $triplet_lib) = add_transcript_to_phase_lib($phase_lib, $triplet_lib, $exon_struct, $strand, $max_tr_rank);
    }
    
    #Disconnect
    $dbh_ens->disconnect;
    
    return ($phase_lib, $triplet_lib);
}

#Add transcript to phase lib
sub add_transcript_to_phase_lib{
    
    #Catch
    my $phase_lib = $_[0];
    my $triplet_lib = $_[1];
    my $exon_struct = $_[2];
    my $strand = $_[3];
    my $max_tr_rank = $_[4];
    
    my $cur_phase = 0;
    my $cur_transcriptomic_pos = 1;
    my $cur_triplet = "";
    if($strand eq '1'){
        for(my $i=1;$i<=$max_tr_rank;$i++){
            my $position = $exon_struct->{$i}->{'start'};
            while($position<=$exon_struct->{$i}->{'end'}){
                $phase_lib->{$strand}->{$position}->{"phase"} = $cur_phase; #Phase lib
                #Get new triplet for phase 0
                if($cur_phase==0){
                    $cur_triplet = substr($exon_struct->{'sequence'},$cur_transcriptomic_pos-1,3);
                }
                $triplet_lib->{$strand}->{$position} = $cur_triplet;#Add triplet to triplet lib
                $phase_lib->{$strand}->{$position}->{"transcriptomic_pos"} = $cur_transcriptomic_pos; #For rel position calculation
                #Get ready for next position
                $cur_phase++;
                if($cur_phase == 3){
                    $cur_phase = 0;
                }
                $cur_transcriptomic_pos++;
                $position++;
            }
        }
        #Save the last position as the sequence length
        my $seq_length = $cur_transcriptomic_pos;
        #Then use this length to calculate normalized positions in the sequence
        for(my $i=1;$i<=$max_tr_rank;$i++){
            my $position2 = $exon_struct->{$i}->{'start'};
            while($position2<=$exon_struct->{$i}->{'end'}){
                $phase_lib->{$strand}->{$position2}->{"transcriptomic_pos"} = $phase_lib->{$strand}->{$position2}->{"transcriptomic_pos"} / $seq_length;
                $position2++;
            }
        }
    } elsif ($strand eq '-1'){
        for(my $i=1; $i<=$max_tr_rank; $i++){
            my $position = $exon_struct->{$i}->{'start'}; #strand relative
            while($position>=$exon_struct->{$i}->{'end'}){
                $phase_lib->{$strand}->{$position}->{"phase"} = $cur_phase;
                if($cur_phase==0){
                    $cur_triplet = substr($exon_struct->{'sequence'},$cur_transcriptomic_pos-1,3);
                }
                $triplet_lib->{$strand}->{$position} = $cur_triplet;
                $phase_lib->{$strand}->{$position}->{"transcriptomic_pos"} = $cur_transcriptomic_pos;
                #Set for next position
                $cur_phase++;
                if($cur_phase == 3){
                    $cur_phase = 0;
                }
                
                $cur_transcriptomic_pos++;
                $position--;
            }
        }
        #Save the last position as the sequence length
        my $seq_length = $cur_transcriptomic_pos;
        #Calculate normalized positions
        for(my $i=1; $i<=$max_tr_rank; $i++){
            my $position = $exon_struct->{$i}->{'start'}; #strand relative
            while($position>=$exon_struct->{$i}->{'end'}){
                $phase_lib->{$strand}->{$position}->{"transcriptomic_pos"} = $phase_lib->{$strand}->{$position}->{"transcriptomic_pos"} / $seq_length;
                $position--;
            }
        }
    }
    
    return $phase_lib, $triplet_lib;
}

## Construct exon structure for certain transcript
sub get_exon_struct_transcript{
    
    #Catch
    my $dbh = $_[0];
    my $transcript = $_[1];
    my $chr = $_[2];
    
    #Init
    my $exon_struct = {};
    
    #Select transcript attributes
    my $query = "SELECT tr.start_exon_id, tr.seq_start, tr.end_exon_id, tr.seq_end, t.seq_region_strand FROM transcript as t JOIN translation as tr ON t.canonical_translation_id=tr.translation_id WHERE t.transcript_id='$transcript';";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    my $transcript_atts = $sth->fetchrow_hashref();
    
    #Get ranked exon structure of transcript
    $query = "SELECT et.rank, e.exon_id, e.seq_region_start, e.seq_region_end, e.phase, e.end_phase FROM exon_transcript as et JOIN exon as e ON et.exon_id=e.exon_id WHERE et.transcript_id='$transcript';";
    $sth = $dbh->prepare($query);
    $sth->execute();
    
    my $ranked_exons = $sth->fetchall_hashref('rank');
    
    #Get max rank
    $query = "SELECT MAX(rank) FROM exon_transcript WHERE transcript_id='$transcript';";
    $sth = $dbh->prepare($query);
    $sth->execute();
    my @result = $sth->fetchrow_array();
    
    my $max_rank = $result[0];
    my $tr_rank = 1;
    my $start_exon_passed = 'False';
    my $max_tr_rank=0;
    
    if($transcript_atts->{'seq_region_strand'}==1){
        for(my $rank=1;$rank<=$max_rank;$rank++){
            #Search for the translation start exon
            if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'start_exon_id'}){
                #Correct for translation start site
                $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_start'} + $transcript_atts->{'seq_start'} - 1;
                #Check if start exon is also the end exon
                if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'end_exon_id'}){
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_start'} + $transcript_atts->{'seq_end'} - 1;
                    $max_tr_rank=$tr_rank;
                } else {
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_end'};
                }
                $exon_struct->{$tr_rank}->{'seq'} = get_sequence($chr,$exon_struct->{$tr_rank}->{'start'},$exon_struct->{$tr_rank}->{'end'});
                #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                $start_exon_passed = 'True';
                $tr_rank++;
            } elsif ($start_exon_passed eq 'True'){
                #Check if in last translated exon
                if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'end_exon_id'}){
                    $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_start'};
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_start'} + $transcript_atts->{'seq_end'} - 1;
                    $exon_struct->{$tr_rank}->{'seq'} = get_sequence($chr,$exon_struct->{$tr_rank}->{'start'},$exon_struct->{$tr_rank}->{'end'});
                    #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                    #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                    $max_tr_rank=$tr_rank;
                    last;
                } else {
                    $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_start'};
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_end'};
                    $exon_struct->{$tr_rank}->{'seq'} = get_sequence($chr,$exon_struct->{$tr_rank}->{'start'},$exon_struct->{$tr_rank}->{'end'});
                    #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                    #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                    $tr_rank++;
                }
            }
        }
        $exon_struct->{'sequence'} = "";
        #Get total translated sequence by concatenation
        for(my $rank = 1; $rank<=$max_tr_rank;$rank++){
            $exon_struct->{'sequence'} = $exon_struct->{'sequence'}.$exon_struct->{$rank}->{'seq'};
        }
        #print "Max tr rank: ".$max_tr_rank."\n";
        #print "Total: ".$exon_struct->{'sequence'}."\n\n";
    } elsif($transcript_atts->{'seq_region_strand'}==-1){
        for(my $rank=1;$rank<=$max_rank;$rank++){
            #Search for the translation start exon
            if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'start_exon_id'}){
                $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_end'} - $transcript_atts->{'seq_start'} + 1;
                #Check if start exon is also the end exon
                if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'end_exon_id'}){
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_end'} - $transcript_atts->{'seq_end'} + 1;
                    $max_tr_rank=$tr_rank;
                } else {
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_start'};
                }
                $exon_struct->{$tr_rank}->{'seq'} = revdnacomp(get_sequence($chr,$exon_struct->{$tr_rank}->{'end'},$exon_struct->{$tr_rank}->{'start'}));
                #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                $start_exon_passed = 'True';
                $tr_rank++;
            } elsif ($start_exon_passed eq 'True'){
                #Check if in last translated exon
                if($ranked_exons->{$rank}->{'exon_id'} eq $transcript_atts->{'end_exon_id'}){
                    $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_end'};
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_end'} - $transcript_atts->{'seq_end'} + 1;
                    $exon_struct->{$tr_rank}->{'seq'} = revdnacomp(get_sequence($chr,$exon_struct->{$tr_rank}->{'end'},$exon_struct->{$tr_rank}->{'start'}));
                    #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                    #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                    $max_tr_rank=$tr_rank;
                    last;
                } else {
                    $exon_struct->{$tr_rank}->{'start'} = $ranked_exons->{$rank}->{'seq_region_end'};
                    $exon_struct->{$tr_rank}->{'end'} = $ranked_exons->{$rank}->{'seq_region_start'};
                    $exon_struct->{$tr_rank}->{'seq'} = revdnacomp(get_sequence($chr,$exon_struct->{$tr_rank}->{'end'},$exon_struct->{$tr_rank}->{'start'}));
                    #print "exon: ".$exon_struct->{$tr_rank}->{'seq'}."\n";
                    #print "\tstart: ".$exon_struct->{$tr_rank}->{'start'}."\tend: ".$exon_struct->{$tr_rank}->{'end'}."\n";
                    $tr_rank++;
                }
            }
        }
        $exon_struct->{'sequence'} = "";
        #Get total translated sequence by concatenation
        for(my $rank = 1; $rank<=$max_tr_rank;$rank++){
            $exon_struct->{'sequence'} = $exon_struct->{'sequence'}.$exon_struct->{$rank}->{'seq'};
        }
        #print "Total: ".$exon_struct->{'sequence'}."\n\n";
    }
    
    my $strand = $transcript_atts->{'seq_region_strand'};
    
    return $exon_struct, $strand, $max_tr_rank;
}


##Run plastid to get p site offsets
sub run_plastid{
    
    #Catch
    my $bam = $_[0];
    my $TMP = $_[1];
    my $version = $_[2];
    my $spec = $_[3];
    my $assembly = $_[4];
    my $exp_name = $_[5];
    my $min_l = $_[6];
    my $max_l = $_[7];
    
    print "PLASTID\n";
    
    #Check bam file. If not mentioned, created bam file out of sam file
    if ($bam eq "convert"){
        my $bam_adress = $TMP."/mappingqc/".$samFileName.".bam";
        if (! -e $bam_adress){
            print "Convert sam file to bam file with samtools view\n";
            system("samtools view -bS -o ".$bam_adress." ".$sam);
        } else {
            print "Bam file for plastid already present\n";
        }
        $bam = $bam_adress;
    }
    
    #Get genes.gtf file
    if (! -e $TMP."/Genes"){
        system("mkdir ".$TMP."/Genes");
    }
    if (! -e $TMP."/Genes/genes.gtf"){
        print "Download genes gtf file for plastid\n";
        if(uc($species) eq "SL1344"){
            system("wget -q ftp://ftp.ensemblgenomes.org/pub/release-".$version."/bacteria/gtf/bacteria_23_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344.".$assembly.".".$version.".gtf.gz");
            system("mv Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344.".$assembly.".".$version.".gtf.gz ".$TMP."/Genes/genesTmp.gtf.gz");
            system("gunzip ".$TMP."/Genes/genesTmp.gtf.gz");
        } else {
            system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/gtf/".lc($spec)."//".$spec.".".$assembly.".".$version.".gtf.gz");
            system("mv ".$spec.".".$assembly.".".$version.".gtf.gz ".$TMP."/Genes/genesTmp.gtf.gz");
            system("gunzip ".$TMP."/Genes/genesTmp.gtf.gz");
        }

        #Strip off first lines
        open (GENESTMP,"<", $TMP."/Genes/genesTmp.gtf") || die "Cannot open genes tmp input\n";
        open (GENES,">", $TMP."/Genes/genes.gtf") || die "Cannot open genes gtf output";
        while(my $line = <GENESTMP>){
            if (!($line =~ m/^.{1}!/)){
                print GENES $line;
            }
        }
        close(GENESTMP);
        close(GENES);
        system("rm -rf ".$TMP."/Genes/genesTmp.gtf");
        
        #For SL1344, simulate UTRs
        if(uc($species) eq "SL1344"){
            system("python ".$tool_dir."/simulate_utr_for_prokaryotes.py ".$TMP."/Genes/genes.gtf > ".$TMP."/Genes/genes_with_UTR.gtf");
            system("mv ".$TMP."/Genes/genes_with_UTR.gtf ".$TMP."/Genes/genes.gtf");
        }
        
    } else {
        print "Genes gtf file for plastid already present\n";
    }
    
    #Generate tmp folder for plastid
    if (! -e $TMP."/plastid"){
        system("mkdir ".$TMP."/plastid");
    }
    
    #Generate metagene
    if (! -e $TMP."/plastid/".$exp_name."_rois.txt"){
        print "Generate metagene\n";
        my $command_metagene = "metagene generate -q ".$exp_name." --landmark cds_start --annotation_files ".$TMP."/Genes/genes.gtf 2> /dev/null";
        system($command_metagene);
        system("mv ".$exp_name."_rois.bed ".$TMP."/plastid");
        system("mv ".$exp_name."_rois.txt ".$TMP."/plastid");
    } else {
        print "Plastid metagene already present\n";
    }
    
    #Calculate offsets
    if (! -e $TMP."/plastid/".$exp_name."_p_offsets.txt"){
        print "Index bam file\n";
        system("samtools index ".$bam);
        
        print "Calculate offsets with plastid\n";
        my $command_offsets = "psite -q ".$TMP."/plastid/".$exp_name."_rois.txt ".$exp_name." --min_length ".$min_l." --max_length ".$max_l." --require_upstream --count_files ".$bam." 2> /dev/null";
        system($command_offsets);
        system("mv ".$exp_name."_metagene_profiles.txt ".$TMP."/plastid");
        system("mv ".$exp_name."_p_offsets.txt ".$TMP."/plastid");
        system("mv ".$exp_name."_p_offsets.png ".$TMP."/plastid");
    } else {
        print "Plastid offsets already present\n";
    }
    
    #Read in offsets
    my $offset_hash = {};
    
    open(OFFSET, "<", $TMP."/plastid/".$exp_name."_p_offsets.txt");
    while(my $line = <OFFSET>){
        if($line =~ /^(\d+)\s+(\d+)$/){
            my $length = $1;
            my $offset = $2;
            $offset_hash->{$length} = $offset;
        }
    }
    close(OFFSET);
    
    #Adapt raw plastid offsets
    my $new_offset_hash = {};
    my ($last) = sort {$a<=>$b} values %{$offset_hash};#To determine offsets which were given value 50 in plastid, start for the first with minimal offset value
    for my $key (sort {$a<=>$b} keys %{$offset_hash}) {
        my $length = $key;
        my $offset = $offset_hash->{$key};
        if ($offset==50){
            $offset = $last #Take offset of last RPF length to correct for values given value 50
        }
        $last = $offset;
        $new_offset_hash->{$length} = $offset; #Save in new offset hash
    }
    $offset_hash = $new_offset_hash; #replace
    
    #Save min and max RPF length of offset hash
    $offset_hash->{'min'} = $min_l;
    $offset_hash->{'max'} = $max_l;

    return $offset_hash;
}

##Write offsets to tmp folder for html file
sub offsets_to_csv {
    
    #Catch
    my $offset_hash = $_[0];
    my $TMP = $_[1];
    
    my $outfile = $TMP."/mappingqc/mappingqc_offsets.csv";
    open(my $fw, '>', $outfile);
    
    for(my $key=$offset_hash->{'min'}; $key<=$offset_hash->{'max'}; $key++){
        my $length = $key;
        my $offset = $offset_hash->{$key};
        my $line = $length.",".$offset."\n";
        print $fw $line;
    }
    
    close($fw);
}

### SPLIT SAM PER CHR ###
sub split_SAM_per_chr {
    
    # Catch
    my %chr_sizes = %{$_[0]};
    my $work_dir = $_[1];
    my $sam = $_[2];
    my $unique = $_[3];
    my $mapper = $_[4];
    
    my @splitsam = split(/\//, $sam );
    my $samFileName = $splitsam[$#splitsam];
    @splitsam = split(/\./,$samFileName);
    $samFileName = $splitsam[0];
    
    my ($chr,@mapping_store,$file_in_loc,$file_in,$file_out);
    
    #Remove eventual existing samfiles
    system("rm -f ".$TMP."/mappingqc/".$samFileName."_*.sam");   # Delete existing
    
    # Touch per chr
    foreach $chr (keys %chr_sizes){
        system("touch ".$TMP."/mappingqc/".$samFileName."_".$chr.".sam");   # Touch new
    }
    
    ## Split files into chromosomes
    
    # Open
    open (I,"<".$sam) || die "Cannot open ".$sam." file\n";
    
    
    #For unsorted SAM file (genomic location)
    my $prev_chr="0";
    my $lines = 0;
    my $count_uniq = 0;
    
    while(my $line=<I>){
        
        #Skip annotation lines
        if ($line =~ m/^@/) { next; }
        #Count the amount of reads to check unique/non-unique mapping
        $lines++;
        
        #Process alignment line
        @mapping_store = split(/\t/,$line);
        $chr = $mapping_store[2];
        
        #For STAR: flag 255 means that there is only 1 alignment
        #For TopHat2/HiSat2: Use NH tag
        if ($unique eq "Y"){
            next unless (($mapping_store[4] == 255 && uc($mapper) eq "STAR") || ($line =~ m/NH:i:1\D/ && (uc($mapper) eq "TOPHAT2" || uc($mapper) eq "HISAT2")));
        } elsif ($unique eq "N"){
            #Check the amount of unique reads to be sure SAM file does not come from unique mapping
            if (($mapping_store[4] == 255 && uc($mapper) eq "STAR") || ($line =~ m/NH:i:1\D/ && (uc($mapper) eq "TOPHAT2" || uc($mapper) eq "HISAT2"))){
                $count_uniq++;
            }
            
            #Keep all mappings, also MultipleMapping locations are available (alternative to pseudogenes mapping) GM:07-10-2013
            #Note that we only retain the up until <16 multiple locations (to avoid including TopHat2 peak @ 16)
            #For STAR the maxMultiMap is not included in the output (e.g. if maxmultimap is set to 16, up untill 15 is included)
            #For TopHat2 the maxMultiMap is included in the output (e.g. if maxmultimap is set to 16, up untill 16 is included)
            #In order to have the same actual limit, the maxmultimap is discarded (see record below)
            next unless ( $line !~ m/NH:i:$maxmultimap/ );
        }
        
        # Write off
        if ($prev_chr ne $chr) {
            if ($prev_chr ne "0") { close(A);}
            $file_out = $TMP."/mappingqc/".$samFileName;
            open (A,">>".$file_out."_".$chr.".sam") || die "Cannot open the sep file";
            print A $line;
        }
        elsif ($prev_chr eq $chr) {
            print A $line;
        }
        $prev_chr = $chr;
        
    }
    
    #Check for unique mapping if non-unique option selected by checking if there as many lines as
    if (($lines == $count_uniq) && $unique eq "N"){
        print "\nWARNING: you selected non-unique mappingQC for a file that probably comes from a unique mapping!!\n\n";
    }
    
    
    
    # Close
    close(A);
    close(I);
}

### Create Bin Chromosomes ##

sub create_BIN_chromosomes {
    
    # Catch
    my $BIN_chrom_dir   =   $_[0];
    my $cores           =   $_[1];
    my $chrs            =   $_[2];
    my $work_dir        =   $_[3];
    my $TMP             =   $_[4];
    
    # Create BIN_CHR directory
    if (! -e $BIN_chrom_dir){
        system ("mkdir -p ".$BIN_chrom_dir);
    }
    
    # Create binary chrom files
    ## Init multi core
    my $pm = new Parallel::ForkManager($cores);
    
    foreach my $chr (keys %{$chrs}){
        
        ## Start parallel process
        $pm->start and next;
        
        if (! -e $BIN_chrom_dir."/".$chr.".fa"){
            open (CHR,"<".$TMP."/Chromosomes/".$chr.".fa") || die "Cannot open chr fasta input\n";
            open (CHR_BIN,">".$TMP."/".$chr.".fa");
            
            while (<CHR>){
                #Skip first header line
                chomp($_);
                if ($_ =~ m/^>/) { next; }
                print CHR_BIN $_;
            }
            
            close(CHR);
            close(CHR_BIN);
            system ("mv ".$TMP."/".$chr.".fa ".$BIN_chrom_dir."/".$chr.".fa");
        } else {
            print "\tBinary chromosome file chromosome ".$chr." is already present\n";
        }
            
        $pm->finish;
    }
    
    # Finish all subprocesses
    $pm->wait_all_children;
    
}

## Get coord system id ##
sub get_coord_system_id{
    # Catch
    my $db_ensembl = $_[0];
    my $assembly = $_[1];
    
    my $user = "";
    my $pw = "";
    
    # Connect to ensembl sqlite database
    my $dbh  = DBI->connect('DBI:SQLite:'.$db_ensembl,$user,$pw,
    { RaiseError => 1},) || die "Database connection not made: $DBI::errstr";
    
    # Get correct coord_system_id
    my $query = "";
    if (uc($species eq "SL1344")){
        $query = "SELECT * FROM coord_system WHERE name='chromosome' and version='".$assembly."';";
    } else {
        $query = "SELECT coord_system_id FROM coord_system WHERE name = 'chromosome' AND version = '$assembly'";
    }
    my $execute = $dbh->prepare($query);
    $execute->execute();
    
    my $coord_system_id;
    while(my @result = $execute->fetchrow_array()){
        $coord_system_id = $result[0];
    }
    
    $execute->finish();
    
    # Disconnect
    $dbh->disconnect();
    
    # Return
    return($coord_system_id);
} # Close sub

### GET CHR SIZES ###
sub get_chr_sizes {
    
    # Catch
    my $chromosome_sizes = $_[0];
    
    # Work
    my %chr_sizes;
    open (Q,"<".$chromosome_sizes) || die "Cannot open chr sizes input\n";
    while (<Q>){
        my @a = split(/\s+/,$_);
        $chr_sizes{$a[0]} = $a[1];
    }
    
    return(\%chr_sizes);
}

### GET CHRs ###

sub get_chrs {
    
    # Catch
    my $db          =   $_[0];
    my $us          =   $_[1];
    my $pw          =   $_[2];
    my $chr_sizes   =   $_[3];
    my $assembly    =   $_[4];
    
    # Init
    my $chrs    =   {};
    my $dbh     =   dbh($db,$us,$pw);
    my ($line,@chr,$coord_system_id,$seq_region_id,@ids,@coord_system);
    
    # Get correct coord_system_id
    my $query = "";
    if(uc($species) eq "SL1344"){
        $query="SELECT * FROM coord_system WHERE name='chromosome' and version='".$assembly."';";
    } else {
        $query = "SELECT coord_system_id FROM coord_system where name = 'chromosome' and version = '".$assembly."'";
    }
    my $sth = $dbh->prepare($query);
    $sth->execute();
    @coord_system = $sth->fetchrow_array;
    $coord_system_id = $coord_system[0];
    $sth->finish();
    
    # Get chrs with seq_region_id
    my $chr;
    #foreach (@chr){
    foreach my $key (keys(%{$chr_sizes})) {
        if (uc($species) eq "FRUITFLY"){
            if($key eq "M"){
                $chr = "dmel_mitochondrion_genome";
            }else{
                $chr = $key;
            }
        } elsif (uc($species) eq "YEAST") {
            if($key eq "MT"){
                $chr = "Mito";
            } else {
                $chr = $key;
            }
        } elsif (uc($species) eq "ZEBRAFISH"){
            if($key eq "MT"){
                $chr = "MtDNA";
            } else {
                $chr = $key;
            }
        } else {
            $chr = $key;
        }
        
        my $query = "SELECT seq_region_id FROM seq_region where coord_system_id = ".$coord_system_id."  and name = '".$chr."' ";
        my $sth = $dbh->prepare($query);
        $sth->execute();
        @ids = $sth->fetchrow_array;
        $seq_region_id = $ids[0];
        $chrs->{$key}{'seq_region_id'} = $seq_region_id;
        $sth->finish();
    }
    
    #Disconnect DBH
    $dbh->disconnect();
    
    # Return
    return($chrs);
    
    
}

## Download one chromosome file
sub downloadChromosomeFasta{
    
    #Catch
    my $chr = $_[0];
    my $spec = $_[1];
    my $version = $_[2];
    my $assembly = $_[3];
    
    #Remove existing
    if ($chr eq "MT" || $chr eq "M"){
        system("rm -rf ".$TMP."/Chromosomes/M.fa");
        system("rm -rf ".$TMP."/Chromosomes/MT.fa");
    } else {
        system("rm -rf ".$TMP."/Chromosomes/".$chr.".fa");
    }
    
    #Download
    if ($chr eq "MT" || $chr eq "M"){
        if(uc($species) eq "FRUITFLY"){
            if ($version>75){
                system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".dna.chromosome.dmel_mitochondrion_genome.fa.gz");
                system("mv ".$spec.".".$assembly.".dna.chromosome.dmel_mitochondrion_genome.fa.gz ".$TMP."/Chromosomes/M.fa.gz");
            } else {
                system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".".$version.".dna.chromosome.dmel_mitochondrion_genome.fa.gz");
                system("mv ".$spec.".".$assembly.".".$version.".dna.chromosome.dmel_mitochondrion_genome.fa.gz ".$TMP."/Chromosomes/M.fa.gz");
            }
            system("gunzip ".$TMP."/Chromosomes/M.fa.gz");
            
        } elsif(uc($species) eq "YEAST"){ #Other name 'Mito' for yeast
            if ($version>75){
                system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".dna.chromosome.Mito.fa.gz");
                system("mv ".$spec.".".$assembly.".dna.chromosome.Mito.fa.gz ".$TMP."/Chromosomes/MT.fa.gz");
            } else {
                system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".".$version.".dna.chromosome.Mito.fa.gz");
                system("mv ".$spec.".".$assembly.".".$version.".dna.chromosome.Mito.fa.gz ".$TMP."/Chromosomes/MT.fa.gz");
            }
            system("gunzip ".$TMP."/Chromosomes/MT.fa.gz");
        } elsif(uc($species) eq "ZEBRAFISH" || uc($species) eq "C.ELEGANS"){ #Other name 'MtDNA' for zebrafish and c.elegans
            if ($version>75){
                system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".dna.chromosome.MtDNA.fa.gz");
                system("mv ".$spec.".".$assembly.".dna.chromosome.MtDNA.fa.gz ".$TMP."/Chromosomes/MT.fa.gz");
            } else {
                system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".".$version.".dna.chromosome.MtDNA.fa.gz");
                system("mv ".$spec.".".$assembly.".".$version.".dna.chromosome.MtDNA.fa.gz ".$TMP."/Chromosomes/MT.fa.gz");
            }
            system("gunzip ".$TMP."/Chromosomes/MT.fa.gz");
        } else {
            if ($version>75){
                system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".dna.chromosome.MT.fa.gz");
                system("mv ".$spec.".".$assembly.".dna.chromosome.MT.fa.gz ".$TMP."/Chromosomes/MT.fa.gz");
            } else {
                system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".".$version.".dna.chromosome.MT.fa.gz");
                system("mv ".$spec.".".$assembly.".".$version.".dna.chromosome.MT.fa.gz ".$TMP."/Chromosomes/MT.fa.gz");
            }
            system("gunzip ".$TMP."/Chromosomes/MT.fa.gz");
        }
    } else {
        if(uc($species) eq "SL1344"){
            system("wget -q ftp://ftp.ensemblgenomes.org/pub/release-".$version."/bacteria/fasta/bacteria_23_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344/dna/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344.".$assembly.".dna.chromosome.".$chr.".fa.gz");
            system("gunzip Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344.".$assembly.".dna.chromosome.".$chr.".fa.gz");
            system("mv Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344.".$assembly.".dna.chromosome.".$chr.".fa ".$TMP."/Chromosomes/".$chr.".fa");
        } else {
            if ($version>75){
                system("wget -q ftp://ftp.ensembl.org/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".dna.chromosome.".$chr.".fa.gz");
                system("mv ".$spec.".".$assembly.".dna.chromosome.".$chr.".fa.gz ".$TMP."/Chromosomes/".$chr.".fa.gz");
            } else {
                system("wget -q ftp://ftp.ensembl.org/ensembl/pub/release-".$version."/fasta/".lc($spec)."/dna//".$spec.".".$assembly.".".$version.".dna.chromosome.".$chr.".fa.gz");
                system("mv ".$spec.".".$assembly.".".$version.".dna.chromosome.".$chr.".fa.gz ".$TMP."/Chromosomes/".$chr.".fa.gz");
            }
            system("gunzip ".$TMP."/Chromosomes/".$chr.".fa.gz");
        }
    }
    print "\t\t*) Chromosome ".$chr." downloaded\n";
    
    return;
}

## Download ChromInfo.txt ##
sub download_chrominfo {

    #Catch
    my $TMP = $_[0];
    my $ucsc = $_[1];

    print "Download chromosomal info file\n";
    
    #Init
    my $chromList = {};
    
    if(uc($species) eq "SL1344"){
        #For bacteria, chromosome sizes can be fetched out of the files where the Ensembl DB is build from
        my $can_version = $version + 53; #Bacteria Ensembl releases are 53 less than the other species
        #first, search coord system id
        system("wget -q ftp://ftp.ensemblgenomes.org/pub/release-".$version."/bacteria/mysql/bacteria_23_collection_core_".$version."_".$can_version."_1/coord_system.txt.gz");
        system("gzip -d coord_system.txt.gz");
        my $coord_system_id = 0;
        open(my $IN_CS, '<',"coord_system.txt");
        while(my $line=<$IN_CS>){
            if ($line =~ m/^(\d+)\t\d+\t\w+\t(\w+)\t1/){
                if($2 eq $assembly){
                    $coord_system_id = quotemeta $1;
                }
            }
        }
        close($IN_CS);
        system("rm -rf coord_system.txt");
        
        #Then, search for the chromosome in seq_region.txt
        system("wget -q ftp://ftp.ensemblgenomes.org/pub/release-".$version."/bacteria/mysql/bacteria_23_collection_core_".$version."_".$can_version."_1/seq_region.txt.gz");
        system("gzip -d seq_region.txt.gz");
        open(my $IN_SR, '<', "seq_region.txt");
        while(my $line=<$IN_SR>){
            if ($line =~ m/^\d+\t(\w+)\t$coord_system_id\t(\d+)/){
                $chromList->{$1} = $2;
            }
        }
        close($IN_SR);
        system("rm -rf seq_region.txt");
    } else {
        system("wget -q ftp://hgdownload.cse.ucsc.edu/goldenPath/".$ucsc."/database/chromInfo.txt.gz");
        system("gzip -d chromInfo.txt.gz");
        system("mv chromInfo.txt ".$TMP);
        system("mv ".$TMP."/chromInfo.txt ".$TMP."/tmpChromInfo.txt");
        open(my $IN,'<',$TMP."/tmpChromInfo.txt") or die "Cannot read ".$TMP."/tmpChromInfo.txt";
        while (my $ line=<$IN>){
            if ($line =~ m/^chr(\w{1,4})\t(\d+)\t\S+\n$/){
                if ($species eq "fruitfly"){
                    $chromList->{$1} = $2; # For fruitfly we want to use 'M' as mitochondrial symbol
                } else {
                    if ($1 eq 'M'){
                        $chromList->{'MT'} = $2;
                    } else {
                        $chromList->{$1} = $2;
                    }
                }
            }
        }
        close($IN);
    }
    
    #Write standard chromosome sizes file
    open(my $OUT, '>', $TMP."/ChromInfo.txt");
    foreach my $key (keys(%{$chromList})){
        print $OUT $key."\t".$chromList->{$key}."\n";
    }
    close($OUT);
    system("rm -rf ".$TMP."/tmpChromInfo.txt");

    return;

}

# Given a Cigar string return a double array
# such that each sub array contiains the [ length_of_operation, Cigar_operation]
sub splitCigar {
    my $cigar_string = shift;
    my @returnable;
    my (@matches) = ($cigar_string =~ /(\d+\w)/g);
    foreach (@matches) {
        my @operation = ($_ =~ /(\d+)(\w)/);
        push @returnable, \@operation;
    }
    return \@returnable;
}

### CALCULATE A-SITE OFFSET ###
sub get_offset {
    
    #Catch
    my $len = $_[0];
    my $offset_hash = $_[1];
    
    #Init
    my $offset;
    
    if($len<$offset_hash->{"min"}){
        $offset = $offset_hash->{$offset_hash->{"min"}};
    } elsif($len>$offset_hash->{"max"}){
        $offset = $offset_hash->{$offset_hash->{"max"}};
    } else {
        $offset = $offset_hash->{$len};
    }
    
    return($offset);
}

### Get sequence from chromosomal sequence fasta files ###
sub get_sequence{
    #Catch
    my $chr = $_[0];
    my $start = $_[1];
    my $end = $_[2];
    
    #Init
    my $seq;
    
    #Open chromosomal sequence fasta file (BINARY)
    open(IN, "< ".$BIN_chrom_dir."/".$chr.".fa");
    binmode(IN);
    
    #translate start and end
    my $offset = $start - 1;
    my $length = $end - $start + 1;
    
    #Define reading start position
    seek(IN, $offset, 0);
    #Read in sequence
    read(IN, $seq, $length);
    
    close(IN);
    
    return $seq;
}

### get reverse complement sequence ###
sub revdnacomp {
    #Catch
    my $dna = $_[0];
    
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

#Get transcripts
sub get_can_transcripts{
    
    #Catch
    my $dbh = $_[0];
    my $seq_region_id = $_[1];
    
    #Init
    my $transcripts = [];
    
    #Take all protein-coding genes and take the canonical transcript ID's from them for the chromosome (seq_region_id)
    my $query = "SELECT t.transcript_id FROM transcript as t JOIN gene as g ON g.canonical_transcript_id=t.transcript_id WHERE g.seq_region_id='$seq_region_id' AND g.biotype='protein_coding';";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    while (my @row = $sth->fetchrow_array()){
        push @{$transcripts}, $row[0];
    }
    
    return $transcripts;
}

## GET seq_region_id ##
sub get_seq_region_id{
    
    #Catch
    my $dbh = $_[0];
    my $chr = $_[1];
    my $coord_system_id = $_[2];
    
    #Init
    my $seq_region_id;
    
    #chr nomenclature ens
    if (uc($species) eq "FRUITFLY"){
        if($chr eq "M"){
            $chr = "dmel_mitochondrion_genome";
        }
    } elsif (uc($species) eq "YEAST") {
        if($chr eq "MT"){
            $chr = "Mito";
        }
    } elsif (uc($species) eq "ZEBRAFISH" or uc($species) eq "C.ELEGANS"){
        if($chr eq "MT"){
            $chr = "MtDNA";
        }
    }
    
    my $query = "SELECT seq_region_id FROM seq_region WHERE coord_system_id = '$coord_system_id' AND name = '$chr';";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    $seq_region_id = $sth->fetch()->[0];
    $sth->finish();
    
    return $seq_region_id;
}

### DBH ###
sub dbh {
    
    # Catch
    my $db  = $_[0];
    my $us	= $_[1];
    my $pw	= $_[2];
    
    # Init DB
    my $dbh = DBI->connect($db,$us,$pw,{ RaiseError => 1 },) || die "Cannot connect: " . $DBI::errstr;
    
    return($dbh);
}

### Help text ###
sub print_help_text {
    
    my $help_string = "\n\nMappingQC (Stand-alone version)

    MappingQC is a tool to easily generate some figures which give a nice overview of the quality of the mapping of ribosome profiling data. More specific, it gives an overview of the P site offset calculation, the gene distribution and the metagenic classification. Furthermore, MappingQC does a thorough analysis of the triplet periodicity and the linked triplet phase (typical for ribosome profiling) in the canonical transcript of your data. Especially, the link between the phase distribution and the RPF length, the relative sequence position and the triplet identity are taken into account.
        
    Input parameters:
    --help                  this helpful screen
    --work_dir              working directory to run the scripts in (default: current working directory)
    --experiment_name       customly chosen experiment name for the mappingQC run (mandatory)
    --samfile               path to the SAM/BAM file that comes out of the mapping script of PROTEOFORMER (mandatory)
    --cores                 the amount of cores to run the script on (integer, default: 5)
    --species               the studied species (mandatory)
    --ens_v                 the version of the Ensembl database you want to use
    --tmp                   temporary folder for storing temporary files of mappingQC (default: work_dir/tmp)
    --unique                whether to use only the unique alignments.
    Possible options: Y, N (default Y)
    --mapper                the mapper you used to generate the SAM file (STAR, TopHat2, HiSat2) (default: STAR)
    --maxmultimap           the maximum amount of multimapped positions used for filtering the reads (default: 16)
    --ens_db                path to the Ensembl SQLite database with annotation info. If you want mappingQC to download the right Ensembl database automatically for you, put in 'get' for this parameter (mandatory)
    --offset                the offset determination method.
                                Possible options:
                                - plastid: calculate the offsets with Plastid (Dunn et al. 2016)
                                - standard: use the standard offsets from the paper of Ingolia et al. (2012) (default option)
                                - from_file: use offsets from an input file
    --plastid_bam           the mapping bam file for Plastid offset generation (default: convert)
    --min_length_plastid    the minimum RPF length for Plastid offset generation (default 22)
    --max_length_plastid    the maximum RPF length for Plastid offset generation (default 34)
    --offset_file           the offsets input file
    --min_length_gd         minimum RPF length used for gene distributions and metagenic classification (default: 26).
    --max_length_gd         maximum RPF length used for gene distributions and metagenic classification (default: 34).
    --outfolder             the folder to store the output files (default: work_dir/mQC_output)
    --tool_dir              folder with necessary additional mappingQC tools. More information below in the dependencies section. (default: search for the default tool directory location in the active conda environment)
    --plotrpftool           the module that will be used for plotting the RPF-phase figure
                                Possible options:
                                - grouped2D: use Seaborn to plot a grouped 2D bar chart (default)
                                - pyplot3D: use mplot3d to plot a 3D bar chart. This tool can suffer sometimes from Escher effects, as it tries to plot a 3D plot with the 2D software of pyplot and matplotlib.
                                - mayavi: use the mayavi package to plot a 3D bar chart. This tool only works on local systems with graphical cards.
    --outhtml               custom name for the output HTML file (default: work_dir/mQC_experiment_name.html)
    --outzip                custom name for output ZIP file (default: work_dir/mQC_experiment_name.zip)
    ";
    
    print $help_string."\n";
    
}
