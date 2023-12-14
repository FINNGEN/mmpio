# MMP::io

Tool to merge many GWAS summary statistics into one file ready to use in [MMP](https://geneviz.aalto.fi/MMP/dashboard/).


## For users

**Download**

1. Make and go to a directory for MMP::io files, for example:

   ```bash
   mkdir ~/mmpio
   cd ~/mmpio
   ```

2. Download mmp from the [release page](https://github.com/FINNGEN/mmpio/releases)

3. Move it to the current mmpio directory you created in (1.):

   ```bash
   mv ~/Downloads/mmp-macos-arm64.tar .
   ```

4. Extract the archive:

   ```bash
   tar -xf mmp-macos-arm64.tar
   ```


**Configure**

Still in the mmpio directory you created in (1.), do:

1. Copy `config.json.sample` to `config.json`

2. Edit `config.json` with name, paths and columns for each of your datasets

3. Specify groups of input files to be used for heterogeneity testing.


**Run**

Your current directory should now be the mmpio directory created in (1.) and it should contain the files `mmp` and `config.json`.

Still in the mmpio directory you created in (1.), do:

```bash
./mmp
```

This outputs a `mmp.tsv` file ready for upload on [MMP](https://geneviz.aalto.fi/MMP/dashboard/).


> [!NOTE]
> **macOS users:** You may get a warning from macOS preventing you from running the binary:
>
> > "mmp" cannot be opened because the developer cannot be verified.
>
> To bypass this you need to:
> 1. Open the mmpio directory created in (1.) with Finder.
> 2. Then *Right-click* on the `mmp` binary file, then *Open*. (double-clicking to open will not work)
> 3. Close the terminal window that opened.
> 4. You should now be able to run the mmp binary with `./mmp`


## For developers

**Requirements**

In order to make changes to the mmp program, you need to have the following:

1. The Go programming language.
   Installation instructions here: https://go.dev/doc/install

**Building**

To build the `mmp` binary:

1. Clone this git repository.
2. Go to the `src` directory within this repository.
3. Run `go build -v` to make the `mmp` binary.
