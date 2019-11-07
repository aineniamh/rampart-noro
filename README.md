# rampart-noro

This pipeline complements [``RAMPART``](https://github.com/artic-network/rampart) and continues downstream analysis to consensus level.

## Installing
Clone this repository:

```
git clone https://github.com/aineniamh/rampart-noro.git
```

Create conda environment and activate it:

```
cd rampart-noro
conda env create -f environment.yml
conda activate rampart-noro
```

## Running

Create run folder:

```
mkdir [run_name]
cd [run_name]
```

Where `[run_name]` is whatever you are calling todays run (as specified in MinKNOW).

Run RAMPART:

```
rampart --protocol ../rampart-noro --basecalledPath ~/MinKNOW/data/reads/[run_name]/pass --annotationOptions barcode_set=[native | rapid | pcr | all]
```

`basecalledPath` should be set to whereever MinKNOW/guppy is going to write its basecalled files.



Open a web browser to view [http://localhost:3000](http://localhost:3000)

