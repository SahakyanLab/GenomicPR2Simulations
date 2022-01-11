#!/bins/

pwd=$(pwd)
save_file="png"

# download species files and analysed meta-data
for species in "01-Prokaryotes" "02-Eukaryotes" "03-Viruses"
do
    echo "Running Download_files.R script..."
    Rscript Download_files.R $pwd $species
    if [ "$species" == "03-Viruses" ]
        then
            Rscript CleanData.R $pwd $species
    fi
done

# obtain meta-data plots for each organism group separately
for species in "01-Prokaryotes" "02-Eukaryotes" "03-Viruses"
do
    echo "Obtaining meta-data plots for" $species "in" $save_file "format..."
    Rscript Analysis.R $pwd $save_file $species
done

# obtain meta-data plots for all species
echo "Obtaining meta-data plots for all species in" $save_file "format..."
cd "./00-All_Species/"
pwd=$(pwd)
Rscript Analysis.R $pwd $save_file
cd ../