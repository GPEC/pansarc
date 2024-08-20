# pansarc
Pan sarcoma fusion gene nanostring analyzer

2019-05-13: PANSARC_RSF no longer works.  Please use PANSARC_RSFv2

installation instruction:
1. install R, R studio
2. install jdk
3. install.packages(c("assertthat","xlsx","reshape2","pander","rChoiceDialogs"))
4. install pansarc_*.tar.gz (latest version)
5. install miktex.org: https://miktex.org/download
6. (on Ubuntu or other Linux), need to install pdflatex
    - Ubuntu: > sudo apt install texlive-latex-base texlive-latex-extra

note: installation for miktex talks a long time!!!

5. unzip PANSARC_RSFv2.zip
6. update probe medians on Nanostring probe list with medians.xlsx
7. use RStudio to open PANSARC_RSF.R
8. change working directory to [GitHub repo home]/PANSARC_RSFv2 (e.g. "/home/samuelc/Documents/workspace/R/pansarc/PANSARC_RSFv2")
8. click "source" in RStudio to run the program.

build instruction (pansarc_*.tar.gz package):
1. change working directory to [GitHub repo home]/pansarc (e.g. "/home/samuelc/Documents/workspace/R/pansarc/pansarc")
2. update version number and date in [GitHub repo home]/pansarc/pansarc/DESCRIPTION (only if there are changes made to source codes)
3. (within R) >devtools::build()

note: testing files located at (MAPcore Microsoft Teams) -> Documents -> IT resources -> GPEC server and apps -> pansarc testing

### installation log:
2024-08-20: installed pansarc_1.1.tar.gz and copied PANSARC_RSFv2.zip to Optiplex XE2 (JBRC room 412)
2019-12-10: installed pansarc_0.1.1.tar.gz and copied PANSARC_RSFv2.zip to Optiplex XE2 (JBRC room 412)
2019-05-13: installed pansarc_0.1.1.tar.gz and copied PANSARC_RSFv2.zip to Julie Ho's laptop
2018-05-14: installed pansarc_0.1.1.tar.gz on Dell PC attached to Applied \Biosystems QuantStudio 6 Flex (room 3405)
