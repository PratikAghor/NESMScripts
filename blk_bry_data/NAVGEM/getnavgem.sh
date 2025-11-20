#!/bin/bash
#SBATCH -Jgetnavgem
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=24G
#SBATCH -t96:00:00
#SBATCH -oReport-%j.out
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
# get navgem data
year=2019

# varArray=("uv-10m" "solrad" "dlwsfc" "ttlpcp" "temp2m" "spchum")
# varArray=("dlwsfc" "ttlpcp" "temp2m" "spchum")
varArray=("dlwsfc")

for i in ${!varArray[@]}; do
  echo "var $i is ${varArray[$i]}"
  var=${varArray[$i]} 
  # file_name_one="navgem2.0_0.176c-sea_${year}_03hr_${var}.D"
  if [[ "$var" == "ttlpcp" ]]; then
    file_name_one="navgem1.4_0.281c-std_${year}_03hr_${var}.D"
  elif [[ "$var" == "solrad" ]]; then
    file_name_one="navgem1.4_0.281c-sec_${year}_03hr_${var}.D"
    file_url_one="http://data.hycom.org/datasets/force/NAVGEM/navgem1.4_0.281c/3hourly/${file_name_one}"
  elif [[ "$var" == "uv-10m" ]]; then
    file_name_one="navgem1.4_0.281c-sec_${year}_03hr_uv-10m.D"
    # file_name_one="navgem2.0_0.176c-sea_${year}_03hr_${var}.D"
    file_url_one="http://data.hycom.org/datasets/force/NAVGEM/navgem1.4_0.281c/3hourly/${file_name_one}"
  else
    file_name_one="navgem1.4_0.281c-sea_${year}_03hr_${var}.D"
  fi

  # if [[ "$var" != "solrad" ]]; then
  #	file_url_one="http://data.hycom.org/datasets/force/NAVGEM/navgem2.0_0.176c/3hourly/${file_name_one}"
  # fi

  file_url_one="http://data.hycom.org/datasets/force/NAVGEM/navgem1.4_0.281c/3hourly/${file_name_one}"

  echo "dowloading $var for ${year}"
  wget --user-agent=" Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36" --no-check-certificate -O ${file_name_one} -c ${file_url_one}
  echo "done"
  echo "#########################"
done
