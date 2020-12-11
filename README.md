# MEDUSA

A Deep Learning based protein flexibility prediction tool.


## Abstract

Information on the protein flexibility is essential to understand crucial molecular mechanisms such as protein stability, interactions with other molecules and protein functions in general. B-factor obtained in the X-ray crystallography experiments is the most common flexibility descriptor available for the majority of the resolved protein structures. Since the gap between the number of the resolved protein structures and available protein sequences is exponentially growing, it is important to provide computational tools for protein flexibility prediction from amino acid sequence. In the current study, we report a Deep Learning based protein flexibility prediction tool MEDUSA (https://www.dsimb.inserm.fr/MEDUSA). MEDUSA uses evolutionary information extracted from protein homologous sequences and amino acid physico-chemical properties as input for a convolutional neural network to assign a flexibility class to each protein sequence position. Trained on a non-redundant dataset of X-ray structures, MEDUSA outperforms the current state-of-the-art flexibility prediction tool PROFbval for a binary classification problem and provides more precise flexibility prediction in three and five classes. MEDUSA is freely available as a web-server providing a clear visualization of the prediction results. Analysis of the MEDUSA output allows a user to identify the potentially highly deformable protein regions and general dynamics properties of the protein.

## Dependencies

### 1. Install Docker

**Linux**

Docker provides a convenient script to install itself (Linux: Ubuntu, Debian|Raspbian, CentOS|RHEL)
```bash
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

```bash
$ git clone https://www.dsimb.inserm.fr/git/gabrielctn/medusa.git
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

### Build and run the Docker container

The container is based on the image [tensorflow/tensorflow:2.3.0](https://hub.docker.com/layers/tensorflow/tensorflow/2.3.0/images/sha256-7bc36fe0ca1a051a808122e87f5438614b371263515df4794abef9a78440af8b?context=explore)

```bash
# Build the container and name (tag) it conveniently: medusa
$ docker build -t medusa .
# Run it. Change paths and names accordingly
$ docker run -it --rm \
    --user "$(id -u):$(id -g)" \ # Launch docker as user's id
    -v /path/to/database:/database:ro \ # Bind mount as read-only the database for HHblits
    -v $(pwd):/project \ # Bind mount the project's folder
    medusa \ # The name of the container we launch
    -i sequence.fasta \ # Fasta file containing the target sequence
    -d uniclust30_2016_09 \ # Name of the database for HHblits (c.f --help for more details)
    -o ./output_dir # Directory which will contain results

# One-liner, for convenience
$ docker run -it --user "$(id -u):$(id -g)" -v /path/to/database:/database:ro -v $(pwd):/project medusa -i sequence.fasta -d uniclust30_2016_09 -o output_dir
```

### Help

```bash
$ docker run medusa
or
$ docker run medusa --help

Usage:
        -i    | --seq         (Required)     Input Fasta sequence (file)
        -d    | --database    (Required)     Name of the database for HHBlits.
                                             Ex: UniRef30_2020_06
                                             If your path is '/path/to/Uniclust/UniRef30_2020_03_a3m.ffdata'
                                             provide: -d UniRef30_2020_03
        -o    | --outdir      (Required)     Path to the output directory.
        -h    | --help        (Optionnal)    Brings up this help
```

### Example

```bash
$ docker run -it --user "$(id -u):$(id -g)" -v /home/cretin/uniclust:/database:ro -v $(pwd):/project medusa -i test.seq.fasta -d uniclust30_2016_09 -o ./outdir
Run HHblits ... done
Run HHfilter ... done
Reformat ... done
Transform alignment to frequence matrix ... done
Create AAindex and one-hot encodings ... done
Translate encodings into vectors ... done
Merge vectors ... done
Run predictions ... done

Results can be found in ./outdir/MEDUSA.jy5J-20201208215144

$ ls ./outdir/MEDUSA.jy5J-2020120
logs  medusa_job_results.tar.gz  prediction
```
