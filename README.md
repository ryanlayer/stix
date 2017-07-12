# STIX

## Build
    git clone https://github.com/ryanlayer/giggle.git
    cd giggle
    make
    cd ..
    wget http://www.sqlite.org/2017/sqlite-amalgamation-3170000.zip
    unzip sqlite-amalgamation-3170000.zip
    git clone https://github.com/ryanlayer/stix.git
    make


## Examples
### DEL

#### 19:12694867-12698924
##### STIX
    STIX_ZERO STIX_ONE STIX_QUANTS STIX_QUANT_DEPTHS
    0         0        0,4,11      0,0,1,1
##### 1KG low-coverge SV callset (Sudmant et al., Nature 2015) ([VCF](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/), [Paper](http://www.nature.com/nature/journal/v526/n7571/full/nature15394.html))
    #CHROM POS       REF ALT   INFO                                        NA12878 NA12890 NA12889
    19     12694867  G   <CN0> CIEND=0,1;CIPOS=0,1;END=12698924            1|1     1|1     1|1
##### SVTYPER high-coverage genotypes
             GT  GQ  SQ     GL         DP RO AO QR QA RS AS ASC RP AP AB
    NA12878  1/1 13  796.49 -82,-3,-2  40 9  30 9  30 0  0  0   9  30 0.77    
    NA12890  1/1 10  609.94 -63,-3,-2  31 7  23 7  23 0  0  0   7  23 0.77
    NA12889  0/1 139 794.47 -80,-1,-15 68 33 34 32 33 20 0  0   12 33 0.51
##### LUMPY/SVTYPER high-coverage calls
    #CHROM POS       REF ALT   INFO                                        NA12878 NA12890 NA12889
    19     12694907  N  <DEL>  CIPOS=-10,69;CIEND=-59,4;END=12698928       1/1     1/1     1/1
High coverage | Low coverage
--------------|-------------
<img src="doc/img/19_12694867-12698924_hi.png" style="width: 4in;"/> | <img src="doc/img/19_12694867-12698924_lo.png" style="width: 4in;"/>


---


#### 5:1022803-1025877
##### STIX
    STIX_ZERO STIX_ONE STIX_QUANTS STIX_QUANT_DEPTHS
    3         0        0,0,0       0,0,0,0
##### 1KG low-coverge SV callset (Sudmant et al., Nature 2015) ([VCF](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/), [Paper](http://www.nature.com/nature/journal/v526/n7571/full/nature15394.html))
    #CHROM POS       REF ALT   INFO                                        NA12878 NA12890 NA12889
    5      1022803   C   <CN0> CIEND=-500,1000;CIPOS=-1000,500;END=1025877 1|1     1|1     1|1
##### SVTYPER high-coverage genotypes
             GT  GQ  SQ   GL          DP  RO  AO QR  QA RS AS ASC RP AP AB
    NA12878  0/0 200 0.00 -0,-18,-61  62  62  0  61  0  27 0  0   34 0  0
    NA12890  0/0 200 0.00 -0,-23,-78  78  78  0  78  0  37 0  0   41 0  0
    NA12889  0/0 200 0.00 -0,-33,-109 110 109 0  109 0  55 0  0   54 0  0
##### LUMPY/SVTYPER high-coverage calls
    N/A
High coverage | Low coverage
--------------|-------------
<img src="doc/img/5_1022803-1025877_hi.png" style="width: 4in;"/> | <img src="doc/img/5_1022803-1025877_lo.png" style="width: 4in;"/>


---


#### 4:113985874-113986369
##### STIX
    STIX_ZERO STIX_ONE STIX_QUANTS STIX_QUANT_DEPTHS
    0         0        4,8,10      0,1,1,1
##### 1KG low-coverge SV callset (Sudmant et al., Nature 2015) ([VCF](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/), [Paper](http://www.nature.com/nature/journal/v526/n7571/full/nature15394.html))
    #CHROM POS        REF ALT   INFO                                        NA12878 NA12890 NA12889
    4      113985874  .   <DEL> END=113986369                               0|0     0|0     0|0
##### SVTYPER high-coverage genotypes
             GT  GQ  SQ      GL         DP RO AO QR QA RS AS ASC RP AP AB
    NA12878  1/1 24  1484.47 -152,-6,-3 74 17 56 17 56 0  0  0   17 56 0.77
    NA12890  0/1 105 446.67  -46,-1,-11 42 22 19 22 19 0  0  0   22 19 0.46
    NA12889  0/1 21  1259.40 -129,-3,-6 71 21 49 21 49 0  0  0   21 49 0.7
##### LUMPY/SVTYPER high-coverage calls
    #CHROM POS        REF ALT   INFO                                        NA12878 NA12890 NA12889
    4      113985951  .   <DEL> CIPOS=-12,2;CIEND=-51,12;END=113986325      1/1     0/1     1/1
High coverage | Low coverage
--------------|-------------
<img src="doc/img/4_113985874-113986369_hi.png" style="width: 4in;"/> | <img src="doc/img/4_113985874-113986369_lo.png" style="width: 4in;"/>



### DUP
### INS

### INV
#### 12:12544868-12546613
##### STIX
##### 1KG low-coverge SV callset (Sudmant et al., Nature 2015) ([VCF](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/), [Paper](http://www.nature.com/nature/journal/v526/n7571/full/nature15394.html))
    #CHROM POS        INFO                                                       NA12878 NA12890 NA12889
    12      12544792  END=12546607;CIEND=-21,21;CIPOS=-21,21;AC=4129;AF=0.824481 1|1     1|1     1|1
##### SVTYPER high-coverage genotypes
             GT  SU PE SR GQ  SQ      GL          DP RO AO QR QA RS AS ASC RP AP AB
    NA12878  0/1 58 58 0  114 1478.17 -150,-2,-14 99 39 60 38 59 20 0  2   18 57 0.61
    NA12890  1/1 54 54 0  25  1356.74 -139,-5,-3  67 15 51 15 51 0  0  0   15 51 0.77
    NA12889  0/1 46 46 0  200 1073.26 -108,-1,-22 94 47 46 47 45 26 0  1   21 44 0.49
##### LUMPY/SVTYPER high-coverage calls
    #CHROM POS        INFO                                   NA12878 NA12890 NA12889
    12     12544868   END=12546613;CIPOS=-2,17;CIEND=-15,19; 0/1     1/1     0/1
High coverage | Low coverage
--------------|-------------
<img src="doc/img/12_12544868-12546613_hi.png" style="width: 4in;"/> | <img src="doc/img/12_12544868-12546613_lo.png" style="width: 4in;"/>

---


#### 12:47290448-47309758
##### STIX
##### 1KG low-coverge SV callset (Sudmant et al., Nature 2015) ([VCF](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/), [Paper](http://www.nature.com/nature/journal/v526/n7571/full/nature15394.html))
    #CHROM POS      INFO                                                        NA12878 NA12890 NA12889
    12     47290470 END=47309756;CIEND=-55,55;CIPOS=-55,55;AC=248;AF=0.0495208  0|0     0|0     1|0
##### SVTYPER high-coverage genotypes
             GT  SU PE SR GQ  SQ     GL          DP  RO  AO QR  QA RS AS ASC RP AP AB
    NA12878  0/0 0  0  0  200 0.00   -0,-39,-130 131 131 0  130 0  64 0  0   66 0  0
    NA12890  0/0 0  0  0  200 0.00   -0,-35,-115 116 116 0  115 0  56 0  0   59 0  0
    NA12889  0/1 79 56 23 146 578.17 -65,-7,-59  120 87  32 86  31 38 0  9   48 22 0.26
##### LUMPY/SVTYPER high-coverage calls
    #CHROM POS      INFO                                   NA12878 NA12890 NA12889
    12     47290448 END=47309758;CIPOS=-3,2;CIEND=-10,5;   0/0     0/0     0/1
High coverage | Low coverage
--------------|-------------
<img src="doc/img/12_47290448-47309758_hi.png" style="width: 4in;"/> | <img src="doc/img/12_47290448-47309758_lo.png" style="width: 4in;"/>

#### 2:89161083-89185670
##### STIX
##### 1KG low-coverge SV callset (Sudmant et al., Nature 2015) ([VCF](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/), [Paper](http://www.nature.com/nature/journal/v526/n7571/full/nature15394.html))
    N/A
##### SVTYPER high-coverage genotypes
             GT  SU PE SR GQ  SQ     GL         DP  RO AO QR QA RS AS ASC RP AP AB 
    NA12878  0/1 0  0  0  4   53.14  -6,-0,-1   4   1  2  1  2  1  0  2   0  0  0.67
    NA12890  0/0 0  0  0  200 0      -0,-28,-93 94  94 0  93 0  42 0  0   51 0  0
    NA12889  0/1 68 44 24 143 383.23 -48,-9,-65 115 89 25 88 24 42 0  2   46 22 0.21
##### LUMPY/SVTYPER high-coverage calls
    #CHROM POS      INFO                                   NA12878 NA12890 NA12889
    2      89161083 END=89185670;CIPOS=0,0;CIEND=-9,0      0/1     0/0     0/1
High coverage | Low coverage
--------------|-------------
<img src="doc/img/2_89161083-89185670_hi.png" style="width: 4in;"/> | <img src="doc/img/2_89161083-89185670_lo.png" style="width: 4in;"/>
