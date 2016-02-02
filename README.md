# viridis
Tools for quantifying gene expression in diatoms.  
    
The natural analysis would follow as:   

1. readsCleaner.py: script to call Trimmomatic and clean the reads.

2. readsMapper.py: script to call STAR and map the reads to the genome.

3.1. readsCounter.py: script to call HTSeq to counts the reads per transcript.  
3.2. samplesCorrelationGrapher.py: script to compute and plot the correlation among samples.  

4.1. cufflinksCaller.py: script to call cufflinks and quantify the mapped reads.  
4.2. pcaGrapher.py: script to compute the PCA on expression values in FPKM.  
4.3. sampleMapper.py: script to discriminate samples in the light-growth space.  
6. logValuesMatrixCreator.py: script to transform absolute (FPKM) into relative expression (log2 fold change). 

