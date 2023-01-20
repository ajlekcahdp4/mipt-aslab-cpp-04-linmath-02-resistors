
base_folder="resources"

red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

current_folder=${2:-./}
passed=true

for file in ${current_folder}/${base_folder}/*.dat; do
    echo -n "Testing ${green}${file}${reset} ... "

    # Check if an argument to executable location has been passed to the program
    $1 --nonverbose < $file > ${current_folder}/$base_folder/temp.tmp

    # Compare inputs
    if $3 ${file}.ans ${current_folder}/${base_folder}/temp.tmp; then
        echo "${green}Passed${reset}"
    else
        echo "${red}Failed${reset}"
        passed=false
    fi
done

if ${passed}
then
    exit 0
else
    # Exit with the best number for an exit code
    exit 666
fi