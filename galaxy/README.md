# Galaxy workflows for PyEuk

The two PyEuk workflows (Figure 1 of the manuscript), authored in Galaxy workflow format 2
(gxformat2). Every computational step except read alignment is a `pyeuk` subcommand.

| file | when to use it |
|---|---|
| `pyeuk-standard.gxwf.yml` | the panel is known: `bwa-mem2` -> `define-windows` -> `call-haplotypes` -> `build-sheet` -> `eukaryotyping` -> `cluster` (distance matrix, cluster sweep, HTML report) |
| `pyeuk-panel-inference.gxwf.yml` | the panel is unknown: `bwa-mem2` against the reference genome -> `derive-panel`, which reconstructs the amplicon panel from read coverage peaks. Its output panel feeds the standard workflow. |

Both are gxformat2. Import through the Galaxy UI (Workflows -> Import -> Upload file) or with
`POST /api/workflows {"from_path": ...}`. Do not round-trip them through `yaml.dump`; edit as text.

The reference input on `define-windows`, `call-haplotypes`, and `derive-panel` is a source
conditional: a FASTA from the history, or a Galaxy-cached genome (the `all_fasta` data table),
matching the `bwa-mem2` idiom.

The companion Galaxy **tools** (`pyeuk_typing`) are proposed to the Tool Shed via
[galaxyproject/tools-iuc#8393](https://github.com/galaxyproject/tools-iuc/pull/8393); the
workflows are also proposed to the IWC via
[galaxyproject/iwc#1356](https://github.com/galaxyproject/iwc/pull/1356). They require the
bioconda `pyeuk` package (0.8.1).
