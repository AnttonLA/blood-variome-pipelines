Pipeline to find the Transcription Factor that is likely affected by a given variant.

It currently has two steps:

1. looks up entries of the ReMAP database (https://remap2022.univ-amu.fr/) to find Transcription Factors that
bind on a given variant.
2. Uses PERFECTOS-APE to find the transcription factors whose binding is likely affected by the variant.


It **requires tabix** to perform the lookup on the ReMap metadata table and find all the relevant studies. 


### Future Additions

It would be nice to incorporate [ffq](https://github.com/pachterlab/ffq) to download the relevant data from ReMap. 
