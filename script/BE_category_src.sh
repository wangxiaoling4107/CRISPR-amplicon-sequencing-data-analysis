#!/bin/sh

script_dir=$PWD"/script"
if [ ! -d $script_dir ];then
        mkdir $script_dir
fi

sed '1d' sample.info.csv|while read line
do
	arr=($line)
	condition=${arr[0]}
	sample="${arr[1]}_${arr[2]}"
	sample=$(echo $sample|sed 's/\./_/g')
	code=${arr[3]}
	edit_win=${arr[4]}   #edit window
	wide_win=${arr[5]}  #wide window
	primer=${arr[6]}
	strand=${arr[7]}
	out_file=$sample"_"$condition"_"$code  ##sample name

	if [[ $sample =~ .*"ABE".* ]];then
		if [[ $strand =~ .*"forward".* ]];then
			desired_edit="AG"
		else
			desired_edit="TC"
		fi
	elif [[ $sample =~ .*"CBE".* ]];then
		if [[ $strand =~ .*"forward".* ]];then
			desired_edit="CT"
		else
			desired_edit="GA"
		fi
	elif [[ $sample =~ .*"Cas".* ]];then #just check randomly
		desired_edit="AG"
	elif [[ $sample =~ .*"Mock".* ]];then #just check randomly
		desired_edit="CT"
	fi
	echo "python /home/xwang3/project/bin/BE_category.py -f $PWD/result/CRISPResso_on_${out_file} --edit_window $edit_win --wide_window $wide_win -d $desired_edit">$script_dir"/s_${out_file}.sh"
done
