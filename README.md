TTT stands for Trivial Tangle Traverser. This tool generates traversals(model sequences) through repetitive genomic tangles that somehow matches coverage and the read alignment.

For help run `./TTT.py --help`

Requires python &ge; 3.9 and dataclasses, pulp, ahocorasick, networkx, statistics, and logging python libraries.

<a href="https://www.biorxiv.org/content/10.64898/2026.03.06.710180v1"> Preprint explaining algorithmic details </a>

### Example usage:
```bash
./TTT.py --graph assembly.gfa --alignment reads.gaf --outdir results_dir --boundary-nodes boundary_nodes.tsv --quality-threshold 20
```

<details>
 <br /> 
<summary><b>Will TTT help with a gap in my scaffold?</b></summary>

 
Generally there are three main reasons for gaps in a scaffold:
* Lack of coverage
  
  TTT searches for the "best" path in the assembly graph that traverses the gap. If there's no path because of the coverage gap &mdash; nothing can be done.
  <p>
  <img width="400" height="400" alt="gap" src="https://github.com/user-attachments/assets/695419c0-a7fb-4728-8454-0433bfe66433" />  
    
  <em> Scaffold &lt;utig4-1497[N100000N:scaffold]&lt;utig4-340 &mdash; nothing can be done </em>
</p>
  
* Long homozygous nodes
  
  Such gaps happen because of the read length being shorter than homozygous nodes. 
  Typical structure looks like a sequence of "bubbles" of similar length, interlaced with long homozygous nodes.
  TTT can be run on such tangles. However, usually if those structures left unresolved in the assembly graph (especially if homozygous nodes are longer than ~100kbp homopolymer-compressed) then there's just no information in the read alignments helping to traverse this region, and thus it will be essentially a random guess.
  <p>
  <img width="400" height="400" alt="diploid_simple_tangle" src="https://github.com/user-attachments/assets/650052c6-5f53-43fa-bdfa-464c8a5d6fdb" />
    
  <em> Scaffolds &lt;utig4-1225&lt;utig4-1224[N5000N:ambig_bubble]&gt;utig4-1511&lt;utig4-1513 and &lt;utig4-1226&lt;utig4-1224[N5000N:ambig_bubble]&gt;utig4-1511&lt;utig4-1512. Because of long homozygous nodes utig4-1224 and utig4-1511 there's just no long reads connecting utig4-1228/utig4-1227 with utig4-1225/utig4-1226 or utig4-1512/utig4-1513. TTT will make a random guess, but so can you </em>
</p>

* Complex repeats
  
  TTT was designed for such cases. However there can be no more than 2 haplotypes in the tangle (so rDNA tangles connecting multiple chromosomes are usually unresolvable).
  Also TTT does not scaffold. Thus, you need to know how to pair incoming and outgoing nodes for two haplotype cases.
  <p>
  <img width="400" height="400" alt="haploid tangle" src="https://github.com/user-attachments/assets/6df23394-811a-49a7-8606-993d2b8f1e89" />  

  <em>Gap caused by repeat array </em>
  </p>

  <p>
  <img width="400" height="400" alt="diploid tangle" src="https://github.com/user-attachments/assets/fc8418dd-1391-4032-ac62-0f9881c9a08c" />
  
  <em>Gap caused by large duplication of homozygous region, present in one of the haplotypes</em>
  </p>

</details>

### Required Arguments:
- `--graph`: Path to the GFA file with the graph structure
- `--alignment`: Path to a file with GraphAligner alignment

Instead of those two options one can use `--verkko-output <verkko output directory>` . In this case internal verkko files for HiFi graph, coverage (ONT) and ONT alignments would be used.

- `--outdir` Output directory

- `--boundary-nodes <boundary_nodes_file>` to locate tangle. 
`boundary_nodes_file` should contain tab separated pairs of incoming and outgoing boundary nodes, one pair by line. Also they should be non-repetive and heterozygous in case of 'diploid' tangles.
Boundary nodes should completely separate the tangle from the rest of the graph &mdash; after their removal there should be no path in remaining graph between tangle nodes and any other non-tangle nodes.

<details>
<summary>Example</summary>
<img width="903" height="895" alt="helo_border" src="https://github.com/user-attachments/assets/693575f7-4bd4-44f0-8774-bc78fbf98224" />

For this tangle decent choice of boundary nodes would be
<br /> 
`utig1-10326     utig1-2575`
<br /> 
`utig1-10327     utig1-2574`
<br /> 
</details>
Currently TTT does not support tangles with more than 2 traversing paths (i.e. most of the rDNA tangles in human-like genomes)

<br />

### Output
TTT outputs two files to the `<outdir>` &mdash; `traversal.multiplicities.csv` with estimated multiplicities of tangle nodes (can be used with Bandage); `traversal.gaf` with the resulting path and, if graph .gfa file contained node sequences &mdash; `traversal.hpc.fasta` with a patch sequence. However, when combined with verkko (since verkko's graph is based on homopolymer-compressed sequences), this patch is also homopolymer compressed. To get non-hpc sequence you'll need to rerun verkko providing `traversal.gaf` with `--path` option &mdash; see <a href="https://github.com/marbl/verkko?tab=readme-ov-file#consensus-for-user-provided-paths"> verkko's manual </a> for details.

### Verkko's final graph coverage fix
In verkko up to (and including ) v2.3.* coverage of the short nodes in tangles in the final graph (assembly.homopolymer-compressed.gfa) is deeply flawed. Currently suggested way is to run TTT.py on the same tangle in hifi-only graph (`2-processGraph/unitig-unrolled-hifi-resolved.gfa` within verkko output directory). Usually this provides better results and does not require realigning ONT reads to graph.

You can find how utig4- nodes match to the utig1- graph in `utig42utig1.gaf` after running `./verkko_coverage_fix/utig4_to_utig1.py <assembly_folder> > utig42utig1.gaf`

Alternatively you can update coverage in final graph running this script
`./verkko_coverage_fix/utig4_coverage_updater.py utig42utig1.gaf <assembly_folder>/assembly.homopolymer-compressed.noseq.gfa <assembly_folder>/2-processGraph/unitig-unrolled-hifi-resolved.ont-coverage.csv > utig4_upt.ont-coverage.csv`

and then pass `utig4_upt.ont-coverage.csv` as `--coverage` in TTT.


### Citation:
 -  Antipov D, Chen Y, Sollitto M, Phillippy AM, Formenti G, Koren S. [Automatic Generation of Model Sequences for Complex Regions in Assembly Graphs](https://www.biorxiv.org/content/10.64898/2026.03.06.710180v1). biorxiv, 2026 


