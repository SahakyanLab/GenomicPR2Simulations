#!/bins/

pwd="$(pwd)/"
save_file="png"

# create folders for storage
mkdir -p ../figures/00-All_Species
mkdir -p ../data/{01-Prokaryotes,02-Eukaryotes,03-Viruses}/{All,PR2_compliance,Raw}
mkdir -p ../figures/{01-Prokaryotes,02-Eukaryotes,03-Viruses}

# download species files and analysed meta-data
###############################################
##### warning - this takes days to finish #####
###############################################
# for species in "01-Prokaryotes" "02-Eukaryotes" "03-Viruses"
# do
#     echo "Running script for $species..."
#     if [ "$species" == "03-Viruses" ]
#         then
#             cd "./$species/"
#             Rscript Download_files.R $pwd $species
#             Rscript CleanData.R $pwd $species
#             cd ../
#         else 
#             cd "./$species/"
#             Rscript Download_files.R $pwd $species
#             cd ../
#     fi
# done

# obtain meta-data plots
cd "./00-All_Species/"
pwd="$(pwd)/"
Rscript Analysis.R $pwd $save_file
cd ../