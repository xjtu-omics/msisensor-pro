### 1. Why MSIsensor-pro?
First, MSIsensor-pro exhibits remarkable advantages in terms of both accuracy and computational cost, compared to the current leading NGS-based MSI scoring methods, especially when processing samples with low sequencing depths or low tumor purities. Then, MSIsensor-pro scores MSI on tumor samples without matched normal requirement, enabling detection of MSI status on patient-derived xenografts/organoids, leukemia, and paraffin-embedded samples. In addition, MSIsensor-pro is able to score MSI using as few as 50 microsatellite sites , indicating its potential to compute MSI status in cancer gene panels.
### 2. What is the cut-off value should I use to determine MSI?  ([#4](https://github.com/xjtu-omics/msisensor-pro/issues/4))
In MSIsensor-pro, we did not provide a cut-off, because the cut-off will be different for different sequencing panel. If you want to use MSIsensor-pro in a new panel for only tumor sample MSI detection, please use the msi status from msisensor or gold standard as true value to obtain a cut-off suit for your panel. By the way, it seems that MSI is a continuous value, and you can just use msi score for you research.

### 3. How many normal or MSI-negative samples should I use to build baseline?
You can use at least 20 normal samples form the same sequencer to build the baselines. If you build baseline with less than 20 test samples, we cannot guarantee the accuracy of MSI calling.

### 4. Could MSIsensor-pro be used for MSI calling with RNA-seq data.

Yes, but we do not have an evaluation of MSI calling with RNA-seq data.