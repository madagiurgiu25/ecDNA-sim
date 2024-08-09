# ecDNA-sim

Simulate ecDNA sequence templates (in .bed format).

This simulator aims to cover a large spectrum of ecDNA structure complexity, defined by the seven topologies in the Giurgiu et al. 2024: (i) Simple circularization, (ii) Simple SV's, (iii) Mixed SV's, (iv) Multi-region, (v) Multi-chromosomal, (vi) Duplications, (vii) Foldbacks.

## Install

```
git clone https://github.com/madagiurgiu25/ecDNA-sim
cd ecdna-sim
python -m pip install -e .
```

## Usage

```
import main as m
import simulate as s

# simulate a multi-fragment rearranged ecDNA, allowing inversion, deletions, foldbacks, originating from multiple-chromosomes
conf = {N:3,
        NEIGHBOR: False,
        SINGLE_CHR: False,
        WEIGHT_CHR: [50,50],
        P_DEL_LEFT: 60,
        P_DEL_RIGHT: 60,
        P_INVERT: 40,
        P_DUP: 0,
        P_FOLDBACK: 40,
       }
m.simulate_simple_mix(s.CHRLEN, s.CHRRANGES)
```

Check more examples under `examples/01_Simulation.ipynb`.

## Input specification

- `N` - number of fragments [1-10]
- `NEIGHBOR` - is the next simulate fragment neighboring the pervious fragment [True/False]
- `SINGLE_CHR` - simulate ecDNA originating from a single chromosome [True/False]
- `WEIGHT_CHR` - weight each chromosomal region 
- `P_DEL_LEFT` - allow small deletions on the left side of the fragment with probability p [0-100]
- `P_DEL_RIGHT` - allow small deletions on the right side of the fragment with probability p [0-100]
- `P_INVERT` - allow inversion of the fragment with probability p [0-100]
- `P_DUP` - allow tandem duplication of the fragment with probability p [0-100]
- `P_FOLDBACK` - allow foldback, i.e. next fragment will overlap with the current one [0-100]


## Output format

```
#chr	start	stop	direction	target	coverage	structure	fragment
chr2	25294254	25594454	-	circ_0	500	circ_0	0
chr2	25367069	25444784	+	circ_0	500	circ_0	1
chr2	25367069	25444784	+	circ_0	500	circ_0	2
chr2	25367069	25444784	+	circ_0	500	circ_0	3
chr12	67786508	68089408	+	circ_0	500	circ_0	4
chr2	10198624	10497524	-	circ_0	500	circ_0	5
chr2	10278364	10387589	-	circ_0	500	circ_0	6
chr2	10310394	10346800	-	circ_0	500	circ_0	7
chr2	10338739	10352420	-	circ_0	500	circ_0	8
chr2	10338739	10352420	-	circ_0	500	circ_0	9
```

### Citation

If you use Decoil for your work please cite our paper:

Madalina Giurgiu, Nadine Wittstruck, Elias Rodriguez-Fos, Rocio Chamorro Gonzalez, Lotte Bruckner, Annabell Krienelke-Szymansky, Konstantin Helmsauer, Anne Hartebrodt, Philipp Euskirchen, Richard P. Koche, Kerstin Haase*, Knut Reinert*, Anton G. Henssen*.**Reconstructing extrachromosomal DNA structural heterogeneity from long-read sequencing data using Decoil**. _Genome Research 2024_, DOI: [https://doi.org/10.1101/gr.279123.124](https://doi.org/10.1101/gr.279123.124)

```
@article{Giurgiu2024ReconstructingDecoil,
    title = {{Reconstructing extrachromosomal DNA structural heterogeneity from long-read sequencing data using Decoil}},
    year = {2024},
    journal = {Genome Research},
    author = {Giurgiu, Madalina and Wittstruck, Nadine and Rodriguez-Fos, Elias and Chamorro Gonzalez, Rocio and Brueckner, Lotte and Krienelke-Szymansky, Annabell and Helmsauer, Konstantin and Hartebrodt, Anne and Euskirchen, Philipp and Koche, Richard P. and Haase, Kerstin and Reinert, Knut and Henssen, Anton G.},
    month = {8},
    pages = {gr.279123.124},
    doi = {10.1101/gr.279123.124},
    issn = {1088-9051}
}
```

Paper repository: [https://github.com/henssen-lab/decoil-paper](https://github.com/henssen-lab/decoil-paper)

### License

The code is distributed under BSD-3-Clause license. See [LICENSE](LICENSE) for details.

### Contact

For any questions do no hesitate to contact us.<br/>
Author: Madalina Giurgiu (madalina.giurgiu@charite.de)

