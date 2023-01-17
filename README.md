# MMP::io

Tool to merge many GWAS summary statistics into one file ready to use in [MMP](https://geneviz.aalto.fi/MMP/dashboard/).


## Usage

**Set the configuration**

1. Copy `config.json.sample` to `config.json`

2. Edit `config.json` with name, paths and columns for each of your datasets

**Run**

After downloading the binary from the [release page](https://github.com/FINNGEN/mmpio/releases), just do:

```sh
tar xf mmp-*.tar
./mmp
```

This outputs a `mmp.tsv` file ready for upload on [MMP](https://geneviz.aalto.fi/MMP/dashboard/).
