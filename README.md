# Freefly: Rapid unbiased fragmentomics analysis of cell-free DNA
Version 1.0.0, May 2023<br />
Author: Kun Sun \(sunkun@szbl.ac.cn\)<br />
<br />
Distributed under the
[GNU General Public License v3.0 \(GPLv3\)](https://www.gnu.org/licenses/gpl-3.0.en.html "GPLv3")
for personal and academic usage only.<br />
For detailed information please refer to the license files under `license` directory.

---

## Installation
`Freefly` contains multiple components, which are mainly written in `C++` for GNU Linux/Unix platforms. After uncompressing
the source package, installation is finished. `Freefly` does not require any specific hardware or OS. The current version
has been tested on CentOS v7.5 with Linux kernel v3.10.0.

`Freefly` depends on the following tools:

- [Ktrim](https://github.com/hellosunking/Ktrim "Ktrim")
- [FLASH](http://ccb.jhu.edu/software/FLASH/ "FLASH")

Pre-compiled executables for these tools are included in this packaged directory (compiled with `g++ v4.8.5` and linked
with `libz v1.2.7`. If you could not run them (which is usually caused by low version of `libc++` or `libz` library),
you may re-compile these programs yourself and replace the ones in the `bin` directory (recommended), then re-compile
the in-house progams for `Freefly` via:
```
user@linux$ make clean && make
```

## Run Freefly
The main program is `freefly` under the same directory as this `README.md` file. You can add its path to your `~/.bashrc`
file under the `PATH` variable to call it from anywhere; or you can run the following command to add it to your current
session:
```
user@linux$ PATH=$PATH:$PWD
```

Call `freefly` without any parameters to see the usage (or use '-h' option):
```
Freefly: Rapid unbiased fragmentomics analysis of cell-free DNA

Usage: freefly [options] -1 <R1.fq[.gz]> -2 <R2.fq[.gz]> -o <sid>

Compulsory parameters:
  -1 R1.fq[.gz]     Specify the path to the files containing read 1
  -2 R2.fq[.gz]     Specify the path to the files containing read 2

                    If you have multiple files for your sample, use ',' to separate them;
                    Gzip-compressed files (with a .gz suffix) are supported


  -o out.prefix     Specify the prefix of the output files

Optional parameters:
  -s size           Minimum read size to be kept after trimming (default: 36; must be larger than 10)
  -k kit            Specify the sequencing kit to use built-in adapters
                    Currently supports 'Illumina' (default), 'Nextera', and 'BGI'
  -m kmer           Kmer length to extract motif (default: 4)
  -Q score          Minimum Phred score to extract motif (default: 30)

  -v                Show software version and exit

Freefly is freely available at https://github.com/hellosunking/Freefly
```

### Example
Your data is generated using Illumina platform with 1 lane, you can run `Freefly` using the following command:

```
user@linux$ freefly -1 /path/to/read1.fq.gz -2 /path/to/read2.fq.gz -o test.sample1
```

### Example 2
Your data is generated using BGI platform with 2 lanes, and you want to analyze 6-mer motifs, then you can run
the analysis using the following command:
```
user@linux$ freefly -k BGI -m 6 -1 /path/to/L1.R1.fq.gz,/path/to/L2.R1.fq.gz -2 /path/to/L1.R2.fq.gz,/path/to/L2.R2.fq.gz -o test.sample2
```

## Testing dataset
As cfDNA data contain patients' sensitive information, so we could not provide any test datasets in this package.
For performance evaluations, we suggest the users use public datasets from the literatures,
e.g., [Liang et al. Clin. Transl. Med. 2020](https://db.cngb.org/search/project/CNP0000680/ "Liang et al."), and
[Bujak et al. PLOS Med. 2020](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA578569 "Bujak et al.").

## Outputs explanation
`Freefly` outputs the size distribution and kmer end motifs (including diversity score) of the cfDNA sample.

---
Please send bug reports to Kun Sun \(sunkun@szbl.ac.cn\).<br />
Freefly is freely available at
[https://github.com/hellosunking/Freefly/](https://github.com/hellosunking/Freefly/ "Freefly @ Github").

