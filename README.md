# viridis
Tools for quantifying gene expression in diatoms.  
    
The natural analysis would follow as:   

cufflinksCaller.py: script to call cufflinks and quantify the mapped reads.  
readsCleaner.py: script to call Trimmomatic and clean the reads.  
readsCounter.py: script to call HTSeq to counts the reads per transcript.  
readsMapper.py: script to call STAR and map the reads to the genome.  
samplesCorrelationGrapher.py: script to compute and plot the correlation among samples.
