    git clone https://github.com/samtools/htslib.git
    cd htslib
    autoheader
    autoconf
    ./configure
    make
    cd ..
    git clone https://github.com/ryanlayer/giggle.git
    cd giggle
    make
    cd ..
    wget http://www.sqlite.org/2014/sqlite-amalgamation-3170000.zip
    unzip sqlite-amalgamation-3170000.zip
    git clone https://github.com/ryanlayer/stix.git
    make
