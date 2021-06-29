# MEDUSA

[![Docker Pulls](https://img.shields.io/docker/pulls/dsimb/medusa.svg)](https://hub.docker.com/r/dsimb/medusa)

(C) Yann Vander Meersche, Gabriel Cretin, Alexandre de Brevern, Jean-Christophe Gelly, Tatiana Galochkina  

**A Deep Learning based protein flexibility prediction tool.**


## Abstract

Information on the protein flexibility is essential to understand crucial molecular mechanisms such as protein stability, interactions with other molecules and protein functions in general. B-factor obtained in the X-ray crystallography experiments is the most common flexibility descriptor available for the majority of the resolved protein structures. Since the gap between the number of the resolved protein structures and available protein sequences is exponentially growing, it is important to provide computational tools for protein flexibility prediction from amino acid sequence. In the current study, we report a Deep Learning based protein flexibility prediction tool MEDUSA (https://www.dsimb.inserm.fr/MEDUSA). MEDUSA uses evolutionary information extracted from protein homologous sequences and amino acid physico-chemical properties as input for a convolutional neural network to assign a flexibility class to each protein sequence position. Trained on a non-redundant dataset of X-ray structures, MEDUSA outperforms the current state-of-the-art flexibility prediction tool PROFbval for a binary classification problem and provides more precise flexibility prediction in three and five classes. MEDUSA is freely available as a web-server providing a clear visualization of the prediction results. Analysis of the MEDUSA output allows a user to identify the potentially highly deformable protein regions and general dynamics properties of the protein.

## Dependencies

### 1. Install Docker

**Linux**

Docker provides a convenient script to install itself (Linux: Ubuntu, Debian|Raspbian, CentOS|RHEL)
```term
$ curl -fsSL https://get.docker.com -o get-docker.sh
$ sudo sh get-docker.sh

# Add yourself to Docker group
$ sudo usermod -aG docker <user>
# This will reload your group assignments,
# avoiding the need to logout and log back in
$ su - $USER
```

**Mac**  

[Get Docker Desktop for Mac](https://docs.docker.com/docker-for-mac/install/)  

**Windows**  

[Get Docker Desktop for Windows](https://docs.docker.com/docker-for-windows/install/)  



### Download this repository  

```term
$ git clone https://github.com/DSIMB/medusa.git
$ cd medusa
```

### Database

HHblits requires a sequence database e.g. uniref30.  
If you don't have a sequence database follow the steps explained on HHblits [extensive user-guide](https://github.com/soedinglab/hh-suite/wiki#hh-suite-databases).  

**TL;DR**:  
1. Download the [latest release](http://wwwuser.gwdg.de/~compbiol/uniclust/current_release/) of database from the HHblits repository into an empty directory using command:  
`wget http://wwwuser.gwdg.de/~compbiol/uniclust/[date]]/UniRef[date]_hhsuite.tar.gz`

2. Extract the files using command:  
`tar xzvf UniRef[date]_hhsuite.tar.gz`


### Download the docker image  

You can download the latest build of Medusa docker image (recommended):  

```
$ docker pull dsimb/medusa
```

or build it yourself from the git repository:  

```
$ docker build -t medusa .
```

### Run the Docker container  
  
**IMPORTANT**  
  
| Option | Adapt to your local paths (no trailing slash '/') |   | Do not modify |
|--------|---------------------------------------------------|---|---------------|
| -v     | /path/to/database                                 | : | /database     |
| -v     | $(pwd)                                            | : | /project      |


```term
# Change paths and names accordingly
PATH_DATABASE=/path/to/database
# Path to MEDUSA git repository
PATH_MEDUSA=$(pwd) 

$ docker run -it --rm \
    # Launch docker as user's id
    --user "$(id -u):$(id -g)" \  
    # Bind mount as read-only the database for HHblits (ex: /home/gabriel/UniRef)
    -v ${PATH_DATABASE}:/database:ro \  
    # Bind mount the project's folder
    -v ${PATH_MEDUSA}:/project \  
    # The name of the container we launch
    medusa \  
    # Fasta file containing the target sequence (path relative to project)
    -i ./data/sequence.fasta \  
    # Name of the database prefix for HHblits (c.f --help for more details)
    # If database file is "uniclust30_2016_09_a3m.ffindex", set option to prefix "uniclust30_2016_09"
    -d uniclust30_2016_09 \  
    # Directory which will contain results (path relative to project)
    -o ./results  
```
One-liner, for convenience
```term
$ docker run -it --user "$(id -u):$(id -g)" -v ${PATH_DATABASE}:/database:ro -v ${PATH_MEDUSA}:/project medusa -i ./data/sequence.fasta -d uniclust30_2016_09 -o ./results
```

### Example for running batch prediction  

If you have a multifasta file and you wish to run MEDUSA protein flexibility prediction for all of your sequences,
you can follow these few steps to split your multifasta file into single fasta files and run the docker container
for each of them individually.  

```term
$ ls
multi.fasta
$ cat multi.fasta
>SEQUENCE_1
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG
LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK
IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL
>SEQUENCE_2
SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI
ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH
$ awk '/^>/{s=++d".fasta"} {print > s}' multi.fasta
1.fasta 2.fasta
$ for seq in ./*.fasta; do docker run -it --user "$(id -u):$(id -g)" -v ${PATH_DATABASE}:/database:ro -v ${PATH_MEDUSA}:/project medusa -i ./$seq -d uniclust30_2016_09 -o ./results -c 0 -o 0; done
```


### Help

```term
$ docker run dsimb/medusa
or
$ docker run dsimb/medusa --help

Usage:
        -i    | --seq         (Required)     Path to input Fasta sequence file. The path is relative to the project folder.
        -d    | --database    (Required)     Name of the database for HHBlits.
                                             Ex: UniRef30_2020_06
                                             If your path is '/path/to/Uniclust/UniRef30_2020_03_a3m.ffdata'
                                             please provide: -d UniRef30_2020_03
        -o    | --outdir      (Required)     Path to the output directory.
        -c    | --cpus        (Optionnal)    Number of CPUs to use (for HHblits). Default is 2. Set to 0 for all.
        -m    | --memory      (Optionnal)    Maximum RAM to use in Gb (for HHblits). Default is 3. Set to 0 for all.
        -h    | --help        (Optionnal)    Brings up this help
```

### Example

```term
$ docker run -it --user "$(id -u):$(id -g)" -v /home/cretin/uniclust:/database:ro -v $(pwd):/project dsimb/medusa -i ./data/sequence.fasta -d uniclust30_2016_09 -o ./results -c 0 -m 0
Run HHblits ... done
Run HHfilter ... done
Reformat ... done
Transform alignment to frequence matrix ... done
Create AAindex and one-hot encodings ... done
Translate encodings into vectors ... done
Merge vectors ... done
Run predictions ... done

Results can be found in ./results/MEDUSA.jy5J-20201208215144

$ ls ./results/MEDUSA.jy5J-2020120
logs  medusa_job_results.tar.gz  prediction
```

### Reference  

Vander Meersche, Y., Cretin, G., de Brevern, A. G., Gelly, J. C., & Galochkina, T. (2021). MEDUSA: Prediction of protein flexibility from sequence. Journal of Molecular Biology, 166882. https://doi.org/10.1016/j.jmb.2021.166882

### Issues  

If you encounter any issue, do not hesitate to open an issue.  
