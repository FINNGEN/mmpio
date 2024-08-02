# MMP::io

Tool to merge many GWAS summary statistics into one file ready to use in [MMP](https://geneviz.aalto.fi/MMP/dashboard/).

Documentation for users: [Download](#download), [Configure](#configure) and [Run](#run).

Want to contribute to the code? Chek the section [For developers](#for-developers).

<img width="764" alt="mmpio overview" src="https://github.com/FINNGEN/mmpio/assets/4919840/f97b7196-398c-40eb-955d-c1a66ec04595">



## For users

### Download

1. Make and go to a directory for MMP::io files, for example:

   ```bash
   mkdir ~/mmpio
   cd ~/mmpio
   ```

2. Download mmpio from the [release page](https://github.com/FINNGEN/mmpio/releases)

3. Move it to the current mmpio directory you created in (1.):

   ```bash
   mv ~/Downloads/mmpio-macos-arm64.tar .
   ```

4. Extract the archive:

   ```bash
   tar -xf mmpio-macos-arm64.tar
   ```


### Configure

Still in the mmpio directory you created in (1.), do:

1. Copy `config.json.sample` to `config.json`

2. Edit `config.json` with name, paths and columns for each of your datasets

3. Specify groups of input files to be used for heterogeneity testing.


### Run

Your current directory should now be the mmpio directory created in (1.) and it should contain the files `mmpio` and `config.json`.

Still in the mmpio directory you created in (1.), do:

```bash
./mmpio
```

This outputs a `mmp.tsv` file ready for upload on [MMP](https://geneviz.aalto.fi/MMP/dashboard/).


> [!NOTE]
> Finngen uses ‘23’ to represent the X chromosome.  Please convert any listing of SNPs with X-pos-ref-alt to be 23-pos-ref-alt before running mmpio.
>
> **macOS users:** You may need an extra step to run the downloaded `mmpio` binary due to macOS security settings.
>
> Here is how to allow running the `mmpio` binary:
> ```sh
> xattr -d com.apple.quarantine /path/to/downloaded/mmpio
> ```


## For developers

### Requirements

In order to make changes to the mmpio program, you need to have the following:

1. The Go programming language.
   Installation instructions here: https://go.dev/doc/install

### Building

To build the `mmpio` binary:

1. Clone this git repository.
2. Go to the `src` directory within this repository.
3. Run `go build -v` to make the `mmpio` binary.
