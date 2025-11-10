# Boyer-Moore-Search-with-Alignment-Comparison-Counters-GRCh38-excerpt
A Boyer–Moore (BM) string search that returns match indices plus the number of alignments tried and character comparisons, using the classic bad-character and good-suffix rules. Makes BM mechanics measurable for teaching—see how preprocessing drives fewer alignments and comparisons vs naïve matching.

**Features**   
	•	BoyerMoore preprocessing (bad-character / good-suffix tables).   
	•	Instrumented BM search: returns ([occurrences], n_alignments, n_comparisons).   
	•	Minimal FASTA reader for test sequences.   

**Demo inputs**    
	•	Pattern (example from class):   
GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG   
	•	Text: chr1.GRCh38.excerpt.fasta (concatenated sequence from FASTA).   

**How it works**    
	•	Precompute BC/GS tables via BoyerMoore(p, alphabet).   
	•	Scan right-to-left at each alignment; on mismatch, shift by max(bad_char, good_suffix).   
	•	After a full match, apply match_skip() shift.      
	_   Count:  _   
	  •	Alignments: how many start positions the pattern was aligned to.   
	  •	Comparisons: total character comparisons done across the whole search.   
