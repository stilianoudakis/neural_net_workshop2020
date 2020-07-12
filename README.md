# neural_net_workshop2020
The final project required by a summer workshop in deep learning in R offer by VCU

## Overview

In this project I will be using a fully connected, sequential, neural network to predict TAD boundary regions using ChIP-seq peak defined transcription factor binding sites.

## Data description

All data is self contained in the `preciseTAD` R package. The developmental version can be found [here](https://github.com/stilianoudakis/preciseTAD)

### Response Type

The response vector in this problem is an indicator variable denoting whether or not a called TAD boundary overlaps with a 50 kb width genomic bins (Y=1) or not (Y=0). Thus, this is a binary classification problem. Random under sampling will be performed prior to modelling to create balanced classes.

### Feature space

I will be considering two types of feature spaces for this project:

   + distance-type: distance in bases from the center of the nearest genomic annotation to the center of the bin
   + overlap count-type: the number of peak overlaps within each genomic bin
   
Normalization procedures will be performed specific to each feature space. For overlap count-type, min/max normalization will be implemented. For distance-type, a log2 normalization will be performed.

## Proposed network architecture

I will be proposing a fully connected, sequential neural network model. My network will include 3 hidden layers, each with 10 fully connected nodes. Likewise, I will be including $L_{1}$ and $L_{2}$ regularization, as well as a dropout rate of 20% to account for overfitting. My neural network will be trained on 70% of the full data, with 30% reserved for testing. Additionally, I will be assigning 20% of my training data for validation.