subs_list=($(ls -d ../sub-??/ |cut -d "/" -f 2 |tail -n 4))
initial_bet_f=0.3
single_sub_ref_img=run-01_T1w
remove_temp_files=no
ss_template_folder=templates/mni_icbm152_nlin_sym_09c
ss_template=mni_icbm152_t1_tal_nlin_sym_09c_brain
qwarp_minpatch=21
qwarp_template_blur=0
qwarp_image_blur=3
epi_slice_order=alt+z
epi_slice_ref=1
epi_to_anat_cost=lpc+ZZ
epi_mni_res=3
epi_fwhm=6
num_cpus=8

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$num_cpus
export OMP_NUM_THREADS=$num_cpus

for sub in "${subs_list[@]}"
do

	################################################################
	################### ANATOMICAL PREPROCESSING ###################
	################################################################

	echo "Anatomical preprocessing: $sub"
	mkdir -p ../derivatives/"$sub"/
	
	imgs_list=($(ls ../"$sub"/anat/*.nii.gz | cut -d "/" -f 4 | cut -d "." -f 1))
	
	# Reorder image list based on single_sub_ref_img
	i=0
	
	for a in "${imgs_list[@]}"
	do
	
	    [[ $a == *"$single_sub_ref_img"* ]] && { break; }
	    (( ++i ))
	    
	done
	
	imgs_list=($(echo ${imgs_list[$i]} ${imgs_list[@]:0:$i} ${imgs_list[@]:(($i + 1))}))

	for img in "${imgs_list[@]}"
	do
	
		echo "$img: crop image FOV"
		robustfov \
		-i ../"$sub"/anat/"$img".nii.gz \
		-r ../derivatives/"$sub"/"$img"_rfov.nii.gz \
		-m ../derivatives/"$sub"/"$img"_rfov.mat \
		&> ../derivatives/"$sub"/"$img"_rfov.log
			
		echo "$img: resample to 1mm isotropic"
		flirt \
		-interp sinc \
		-sincwindow hanning \
		-in ../derivatives/"$sub"/"$img"_rfov.nii.gz \
		-ref ../derivatives/"$sub"/"$img"_rfov.nii.gz \
		-applyisoxfm 1 \
		-init ident.mat \
		-out ../derivatives/"$sub"/"$img"_res.nii.gz \
		&> ../derivatives/"$sub"/"$img"_res.log

		echo "$img: initial brain extraction"
		bet \
		../derivatives/"$sub"/"$img"_res.nii.gz \
		../derivatives/"$sub"/"$img"_bet.nii.gz \
		-f "$initial_bet_f" \
		&> ../derivatives/"$sub"/"$img"_bet.log
		
		if [[ "$img" == *"$single_sub_ref_img"* ]]
		then
		
 			echo "$img: ACPC alignment"
			flirt \
  			-dof 12 \
			-in ../derivatives/"$sub"/"$img"_bet.nii.gz \
			-ref templates/mni_icbm152_nlin_sym_09c/mni_icbm152_t1_tal_nlin_sym_09c_brain.nii.gz \
			-omat ../derivatives/"$sub"/"$img"_2acpc.mat \
			&>> ../derivatives/"$sub"/"$img"_2acpc.log
			
			aff2rigid ../derivatives/"$sub"/"$img"_2acpc.mat \
			../derivatives/"$sub"/"$img"_2acpc.mat \
			&>> ../derivatives/"$sub"/"$img"_2acpc.log
			
			flirt \
			-interp sinc \
			-sincwindow hanning \
			-in ../derivatives/"$sub"/"$img"_bet.nii.gz \
			-ref templates/mni_icbm152_nlin_sym_09c/mni_icbm152_t1_tal_nlin_sym_09c_brain.nii.gz \
			-applyxfm \
			-init ../derivatives/"$sub"/"$img"_2acpc.mat \
			-out ../derivatives/"$sub"/"$img"_bet_2acpc.nii.gz \
			&>> ../derivatives/"$sub"/"$img"_2acpc.log

			flirt \
			-interp sinc \
			-sincwindow hanning \
			-in ../derivatives/"$sub"/"$img"_res.nii.gz \
			-ref templates/mni_icbm152_nlin_sym_09c/mni_icbm152_t1_tal_nlin_sym_09c_brain.nii.gz \
			-applyxfm \
			-init ../derivatives/"$sub"/"$img"_2acpc.mat \
			-out ../derivatives/"$sub"/"$img"_2acpc.nii.gz \
			&>> ../derivatives/"$sub"/"$img"_2acpc.log
			
			echo "$img: brain extraction"
			antsBrainExtraction.sh \
			-d 3 \
			-a ../derivatives/"$sub"/"$img"_2acpc.nii.gz \
			-e templates/oasis/T_template0.nii.gz \
			-m templates/oasis/T_template0_BrainCerebellumProbabilityMask.nii.gz \
			-f templates/oasis/T_template0_BrainCerebellumRegistrationMask.nii.gz \
			-o ../derivatives/"$sub"/"$img"_ \
			&> ../derivatives/"$sub"/"$img"_brain_extraction.log
						
			mv ../derivatives/"$sub"/"$img"_BrainExtractionMask.nii.gz \
			../derivatives/"$sub"/"$img"_2acpc_brain_mask.nii.gz
			
			3dcalc -a ../derivatives/"$sub"/"$img"_2acpc.nii.gz \
			-b ../derivatives/"$sub"/"$img"_2acpc_brain_mask.nii.gz \
			-expr 'a*b' \
			-prefix ../derivatives/"$sub"/"$img"_2acpc_brain.nii.gz \
			&> ../derivatives/"$sub"/"$img"_2acpc_brain.log
			
			echo "$img: bias field correction and segmentation"
			antsAtroposN4_alt.sh \
			-d 3 \
			-a ../derivatives/"$sub"/"$img"_2acpc_brain.nii.gz \
			-x ../derivatives/"$sub"/"$img"_2acpc_brain_mask.nii.gz \
			-c 3 \
			-y 2 \
			-y 3 \
			-w 0.25 \
			-o ../derivatives/"$sub"/"$img"_2acpc_brain_norm \
			&> ../derivatives/"$sub"/"$img"_2acpc_brain_norm.log
			
			mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentation0N4.nii.gz \
			../derivatives/"$sub"/"$img"_2acpc_brain_norm.nii.gz
			
			mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentation0N4.nii.gz_BIAS.nii.gz \
			../derivatives/"$sub"/"$img"_2acpc_brain_norm_bias.nii.gz
			
			mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentationPosteriors1.nii.gz \
			../derivatives/"$sub"/"$img"_2acpc_brain_csf_prob.nii.gz
			
			mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentationPosteriors2.nii.gz \
			../derivatives/"$sub"/"$img"_2acpc_brain_gm_prob.nii.gz
			
			mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentationPosteriors3.nii.gz \
			../derivatives/"$sub"/"$img"_2acpc_brain_wm_prob.nii.gz
			
			mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentation.nii.gz \
			../derivatives/"$sub"/"$img"_2acpc_brain_binseg.nii.gz
			
			3dcalc -a ../derivatives/"$sub"/"$img"_2acpc.nii.gz \
			-b ../derivatives/"$sub"/"$img"_2acpc_brain_norm_bias.nii.gz \
			-expr 'a/b' \
			-prefix ../derivatives/"$sub"/"$img"_2acpc_norm.nii.gz \
			&>> ../derivatives/"$sub"/"$img"_2acpc_brain_norm.log
			
			echo "$img: refit files"
			3drefit \
			-view orig \
			-space ORIG \
			../derivatives/"$sub"/"$img"_2acpc_norm.nii.gz \
			&> /dev/null
			
			3drefit \
			-view orig \
			-space ORIG \
			../derivatives/"$sub"/"$img"_2acpc.nii.gz \
			&> /dev/null
			
			3drefit \
			-view orig \
			-space ORIG \
			../derivatives/"$sub"/"$img"_2acpc_brain.nii.gz \
			&> /dev/null			
						
		fi
		
		if [[ "$img" != *"$single_sub_ref_img"* ]]
		then
		
			echo "$img: coregister to ACPC anatomical reference"
			flirt -dof 6 \
			-interp sinc \
			-sincwindow hanning \
			-in ../derivatives/"$sub"/"$img"_bet.nii.gz \
			-ref ../derivatives/"$sub"/"$sub"_"$single_sub_ref_img"_2acpc_brain_norm.nii.gz \
			-out  ../derivatives/"$sub"/"$img"_bet_2acpc.nii.gz \
			-omat ../derivatives/"$sub"/"$img"_bet_2acpc.mat \
			&>> ../derivatives/"$sub"/"$img"_bet_2acpc.log
			
			flirt \
			-interp sinc \
			-sincwindow hanning \
			-in ../derivatives/"$sub"/"$img"_res.nii.gz \
			-ref ../derivatives/"$sub"/"$sub"_"$single_sub_ref_img"_2acpc_norm.nii.gz \
			-applyxfm \
			-init ../derivatives/"$sub"/"$img"_bet_2acpc.mat \
			-out ../derivatives/"$sub"/"$img"_2acpc.nii.gz \
			&>> ../derivatives/"$sub"/"$img"_2acpc.log
			
			echo "$img: brain extraction"
			3dcalc -a ../derivatives/"$sub"/"$img"_2acpc.nii.gz \
			-b ../derivatives/"$sub"/"$sub"_"$single_sub_ref_img"_2acpc_brain_mask.nii.gz \
			-expr 'a*b' \
			-prefix ../derivatives/"$sub"/"$img"_2acpc_brain.nii.gz \
			&> ../derivatives/"$sub"/"$img"_2acpc_brain.log
			
			if [[ "$img" == *"T1w"* ]]
			then
			
				echo "$img: bias field correction and segmentation"
				antsAtroposN4_alt.sh \
				-d 3 \
				-a ../derivatives/"$sub"/"$img"_2acpc_brain.nii.gz \
				-x ../derivatives/"$sub"/"$img"_2acpc_brain_mask.nii.gz \
				-c 3 \
				-y 2 \
				-y 3 \
				-w 0.25 \
				-o ../derivatives/"$sub"/"$img"_2acpc_brain_norm \
				&> ../derivatives/"$sub"/"$img"_2acpc_brain_norm.log
				
				mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentation0N4.nii.gz \
				../derivatives/"$sub"/"$img"_2acpc_brain_norm.nii.gz
			
				mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentation0N4.nii.gz_BIAS.nii.gz \
				../derivatives/"$sub"/"$img"_2acpc_brain_norm_bias.nii.gz
			
				mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentationPosteriors1.nii.gz \
				../derivatives/"$sub"/"$img"_2acpc_brain_csf_prob.nii.gz
			
				mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentationPosteriors2.nii.gz \
				../derivatives/"$sub"/"$img"_2acpc_brain_gm_prob.nii.gz
			
				mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentationPosteriors3.nii.gz \
				../derivatives/"$sub"/"$img"_2acpc_brain_wm_prob.nii.gz
			
				mv ../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentation.nii.gz \
				../derivatives/"$sub"/"$img"_2acpc_brain_binseg.nii.gz
			
				3dcalc -a ../derivatives/"$sub"/"$img"_2acpc.nii.gz \
				-b ../derivatives/"$sub"/"$img"_2acpc_brain_norm_bias.nii.gz \
				-expr 'a/b' \
				-prefix ../derivatives/"$sub"/"$img"_2acpc_norm.nii.gz \
				&>> ../derivatives/"$sub"/"$img"_2acpc_brain_norm.log
			
				echo "$img: refit files"
				3drefit \
				-view orig \
				-space ORIG \
				../derivatives/"$sub"/"$img"_2acpc_norm.nii.gz \
				&> /dev/null									
			
			fi
		
		fi						
		
		if [[ "$remove_temp_files" == *"yes"* ]]
		then
		
			echo "$img: removing temporary files"
			rm -rf ../derivatives/"$sub"/"$img"_bet.nii.gz \
			../derivatives/"$sub"/"$img"_bet_2acpc.nii.gz \
			../derivatives/"$sub"/"$img"_res.nii.gz \
			../derivatives/"$sub"/"$img"_rfov.nii.gz \			
			../derivatives/"$sub"/"$img"_BrainExtractionPrior0GenericAffine.mat \
			../derivatives/"$sub"/"$img"_BrainExtractionBrain.nii.gz \
			../derivatives/"$sub"/"$img"_2acpc_brain_normSegmentationConvergence.txt
			
		fi			

	done
	
	echo "$sub: create average T1w image"
	3dMean -prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
	../derivatives/"$sub"/"$sub"_run-??_T1w_2acpc_norm.nii.gz \
	&> ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.log
	
	echo "$sub: create average brain-extracted T1w image"
	3dMean -prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_brain_norm.nii.gz \
	../derivatives/"$sub"/"$sub"_run-??_T1w_2acpc_brain_norm.nii.gz \
	&> ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_brain_norm.log	
	
	echo "$sub: estimate transformation to standard space"
	3dQwarp -allineate -blur "$qwarp_template_blur" "$qwarp_image_blur" \
	-base "$ss_template_folder"/"$ss_template".nii.gz \
	-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_brain_norm.nii.gz \
	-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm.nii.gz \
	-iwarp \
	-minpatch "$qwarp_minpatch" \
	&> ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm.log

	echo "$sub: apply transformation to standard space"	
	3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
	-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
	-master "$ss_template_folder"/"$ss_template".nii.gz \
	-interp wsinc5 \
	-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_norm.nii.gz \
	&> ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_norm.log
	
	for img in "${imgs_list[@]}"
	do
	
		echo "$img: apply transformation to standard space"	
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$img"_2acpc_brain.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp wsinc5 \
		-prefix ../derivatives/"$sub"/"$img"_2mni_brain.nii.gz \
		&> ../derivatives/"$sub"/"$img"_2mni_brain.log
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$img"_2acpc.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp wsinc5 \
		-prefix ../derivatives/"$sub"/"$img"_2mni.nii.gz \
		&> ../derivatives/"$sub"/"$img"_2mni.log		
	
	done
	
	if [ -f ../derivatives/"$sub"/"$sub"_run-01_T2w_2acpc.nii.gz ]
	then

		echo "$sub: run freesurfer"
		recon-all -i ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-subjid freesurfer \
		-sd ../derivatives/"$sub"/ \
		-T2 ../derivatives/"$sub"/"$sub"_run-01_T2w_2acpc.nii.gz \
		-T2pial \
		-qcache \
		-all &> ../derivatives/"$sub"/"$sub"_freesurfer.log
		
		echo "$sub: convert freesurfer files to nifti/gifti"
		@SUMA_Make_Spec_FS -NIFTI \
		-sid freesurfer \
		-fspath ../derivatives/"$sub"/freesurfer/ \
		&> /dev/null
		
		3dAllineate \
		-input ../derivatives/"$sub"/freesurfer/SUMA/aparc+aseg.nii.gz \
		-master ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_aparc+aseg.nii.gz \
		-final NN \
		-1Dparam_apply '1D: 12@0'\' \
		&> /dev/null
		
		3dAllineate \
		-input ../derivatives/"$sub"/freesurfer/SUMA/aparc+aseg_REN_vent.nii.gz \
		-master ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-prefix ../derivatives/"$sub"/temp_"$sub"_avg_T1w_2acpc_ventricles.nii.gz \
		-final NN \
		-1Dparam_apply '1D: 12@0'\' \
		&> /dev/null
		
		3dcalc -a ../derivatives/"$sub"/temp_"$sub"_avg_T1w_2acpc_ventricles.nii.gz \
		-expr 'step(a)' \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_ventricles.nii.gz \
		&> /dev/null
		
		rm -rf ../derivatives/"$sub"/temp_"$sub"_avg_T1w_2acpc_ventricles.nii.gz
		
		3dcalc -a ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_aparc+aseg.nii.gz \
		-expr 'equals(a,2)+equals(a,41)+equals(a,255)+equals(a,254)+equals(a,253)+equals(a,252)+equals(a,251)' \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_white_matter.nii.gz \
		&> /dev/null
		
		3dAllineate \
		-input ../derivatives/"$sub"/freesurfer/SUMA/lh.ribbon.nii.gz \
		-master ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_lh_cortex.nii.gz \
		-final NN \
		-1Dparam_apply '1D: 12@0'\' \
		&> /dev/null
		
		3dAllineate \
		-input ../derivatives/"$sub"/freesurfer/SUMA/rh.ribbon.nii.gz \
		-master ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_rh_cortex.nii.gz \
		-final NN \
		-1Dparam_apply '1D: 12@0'\' \
		&> /dev/null

		3dcalc -a ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_lh_cortex.nii.gz \
		-b ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_rh_cortex.nii.gz \
		-expr 'step(a+b)' \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_cortex.nii.gz \
		&> /dev/null
		
		3dcalc -a ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_aparc+aseg.nii.gz \
		-expr 'step(equals(a,10)+equals(a,11)+equals(a,12)+equals(a,13)+equals(a,17)+equals(a,18)+equals(a,26)+equals(a,49)+equals(a,50)+equals(a,51)+equals(a,52)+equals(a,53)+equals(a,54)+equals(a,58))' \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_deep_gray_matter.nii.gz \
		&> /dev/null
		
		echo "$sub: transform freesurfer rois to standard space"
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_aparc+aseg.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_aparc+aseg.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_ventricles.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_ventricles.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_white_matter.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_white_matter.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_lh_cortex.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_lh_cortex.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_rh_cortex.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_rh_cortex.nii.gz \
		&> /dev/null

		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_cortex.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_cortex.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_deep_gray_matter.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_deep_gray_matter.nii.gz \
		&> /dev/null
	
	else
	
		echo "$sub: run freesurfer"
		recon-all -i ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-subjid freesurfer \
		-sd ../derivatives/"$sub"/ \
		-qcache \
		-all &> ../derivatives/"$sub"/"$sub"_freesurfer.log
		
		echo "$sub: convert freesurfer files to nifti/gifti"
		@SUMA_Make_Spec_FS -NIFTI \
		-sid freesurfer \
		-fspath ../derivatives/"$sub"/freesurfer/ \
		&> /dev/null
		
		3dAllineate \
		-input ../derivatives/"$sub"/freesurfer/SUMA/aparc+aseg.nii.gz \
		-master ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_aparc+aseg.nii.gz \
		-final NN \
		-1Dparam_apply '1D: 12@0'\' \
		&> /dev/null
		
		3dAllineate \
		-input ../derivatives/"$sub"/freesurfer/SUMA/aparc+aseg_REN_vent.nii.gz \
		-master ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-prefix ../derivatives/"$sub"/temp_"$sub"_avg_T1w_2acpc_ventricles.nii.gz \
		-final NN \
		-1Dparam_apply '1D: 12@0'\' \
		&> /dev/null
		
		3dcalc -a ../derivatives/"$sub"/temp_"$sub"_avg_T1w_2acpc_ventricles.nii.gz \
		-expr 'step(a)' \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_ventricles.nii.gz \
		&> /dev/null
		
		rm -rf ../derivatives/"$sub"/temp_"$sub"_avg_T1w_2acpc_ventricles.nii.gz
		
		3dcalc -a ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_aparc+aseg.nii.gz \
		-expr 'equals(a,2)+equals(a,41)+equals(a,255)+equals(a,254)+equals(a,253)+equals(a,252)+equals(a,251)' \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_white_matter.nii.gz \
		&> /dev/null
		
		3dAllineate \
		-input ../derivatives/"$sub"/freesurfer/SUMA/lh.ribbon.nii.gz \
		-master ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_lh_cortex.nii.gz \
		-final NN \
		-1Dparam_apply '1D: 12@0'\' \
		&> /dev/null
		
		3dAllineate \
		-input ../derivatives/"$sub"/freesurfer/SUMA/rh.ribbon.nii.gz \
		-master ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_norm.nii.gz \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_rh_cortex.nii.gz \
		-final NN \
		-1Dparam_apply '1D: 12@0'\' \
		&> /dev/null

		3dcalc -a ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_lh_cortex.nii.gz \
		-b ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_rh_cortex.nii.gz \
		-expr 'step(a+b)' \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_cortex.nii.gz \
		&> /dev/null
		
		3dcalc -a ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_aparc+aseg.nii.gz \
		-expr 'step(equals(a,10)+equals(a,11)+equals(a,12)+equals(a,13)+equals(a,17)+equals(a,18)+equals(a,26)+equals(a,49)+equals(a,50)+equals(a,51)+equals(a,52)+equals(a,53)+equals(a,54)+equals(a,58))' \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_deep_gray_matter.nii.gz \
		&> /dev/null
		
		echo "$sub: transform freesurfer rois to standard space"
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_aparc+aseg.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_aparc+aseg.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_ventricles.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_ventricles.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_white_matter.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_white_matter.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_lh_cortex.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_lh_cortex.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_rh_cortex.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_rh_cortex.nii.gz \
		&> /dev/null

		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_cortex.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_cortex.nii.gz \
		&> /dev/null
		
		3dNwarpApply -nwarp "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz" \
		-source ../derivatives/"$sub"/"$sub"_avg_T1w_2acpc_deep_gray_matter.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-interp NN \
		-prefix ../derivatives/"$sub"/"$sub"_avg_T1w_2mni_deep_gray_matter.nii.gz \
		&> /dev/null
	
	fi
	
	################################################################
	################### FUNCTIONAL PREPROCESSING ###################
	################################################################
	
	imgs_list=($(ls ../"$sub"/func/*.nii.gz | cut -d "/" -f 4 | cut -d "." -f 1))
	echo "Functional preprocessing: $sub"
	
	for img in "${imgs_list[@]}"
	do
	
		n_tr=$(3dinfo -nt ../"$sub"/func/"$img".nii.gz)
		
		epi_volume_ref=$(echo "$n_tr / 2" |bc)
	
		echo "$img: slice timing correction"
		3dTshift \
		-Fourier \
		-tzero "$epi_slice_ref" \
		-tpattern "$epi_slice_order" \
		-prefix ../derivatives/"$sub"/"$img"-t.nii.gz \
		../"$sub"/func/"$img".nii.gz \
		&> ../derivatives/"$sub"/"$img"-t.log
		
		echo "$img: brain extraction"
		bet ../derivatives/"$sub"/"$img"-t.nii.gz \
		../derivatives/"$sub"/"$img"-tb \
		-F \
		&> ../derivatives/"$sub"/"$img"-tb.log
		
		echo "$img: head motion correction"
		3dvolreg \
		-Fourier \
		-base "$epi_volume_ref" \
		-1Dfile ../derivatives/"$sub"/"$img"-mocopar.1D \
		-1Dmatrix_save ../derivatives/"$sub"/"$img"-mocomatrix.aff12.1D \
		-maxdisp1D ../derivatives/"$sub"/"$img"-mocomaxdisp.1D \
		-prefix ../derivatives/"$sub"/"$img"-tbv.nii.gz \
		../derivatives/"$sub"/"$img"-tb.nii.gz \
		&> ../derivatives/"$sub"/"$img"-tbv.log
		
		3dTstat \
		-prefix ../derivatives/"$sub"/"$img"-tbv_avg.nii.gz \
		../derivatives/"$sub"/"$img"-tbv.nii.gz \
		&> /dev/null
		
		echo "$img: linear epi to anat registration"
		align_epi_anat.py \
		-anat ../derivatives/"$sub"/"$sub"_"$single_sub_ref_img"_2acpc_brain_norm.nii.gz \
		-epi ../derivatives/"$sub"/"$img"-tbv_avg.nii.gz \
		-epi_base 0 \
		-epi2anat \
		-giant_move \
		-anat_has_skull no \
		-epi_strip None \
		-volreg off \
		-tshift off \
		-deoblique off \
		-cost "$epi_to_anat_cost" \
		-output_dir ../derivatives/"$sub"/ \
		&> ../derivatives/"$sub"/"$img"-lin_epi2anat.log
		
		3dAllineate \
		-input ../derivatives/"$sub"/"$img"-tbv_avg_al+orig.HEAD \
		-master ../derivatives/"$sub"/"$sub"_"$single_sub_ref_img"_2acpc_brain_norm.nii.gz \
		-prefix ../derivatives/"$sub"/"$img"-tbv_avg_2acpc.nii.gz \
		-final wsinc5 \
		-1Dparam_apply '1D: 12@0'\' \
		&> /dev/null		
		
		mv ../derivatives/"$sub"/"$img"-tbv_avg_al_mat.aff12.1D \
		../derivatives/"$sub"/"$img"-tbv_avg_2acpc.aff12.1D
		
		rm -rf ../derivatives/"$sub"/"$img"-tbv_avg_al+orig.* \
		../derivatives/"$sub"/"$sub"_"$single_sub_ref_img"_2acpc_brain_norm_al_mat.aff12.1D

		echo "$img: nonlinear epi to anat registration"
		3dQwarp -source ../derivatives/"$sub"/"$sub"_"$single_sub_ref_img"_2acpc_brain_norm.nii.gz \
		-base ../derivatives/"$sub"/"$img"-tbv_avg_2acpc.nii.gz \
		-prefix ../derivatives/"$sub"/"$img"-tbvu_avg_2acpc \
		-lpc -maxlev 0 -verb -iwarp -blur 0 3 \
		&> ../derivatives/"$sub"/"$img"-nonlin_epi2anat.log

		rm -rf ../derivatives/"$sub"/"$img"-tbvu_avg_2acpc_WARP+orig.* \
		../derivatives/"$sub"/"$img"-tbvu_avg_2acpc+orig.*
		
		3dcopy  ../derivatives/"$sub"/"$img"-tbvu_avg_2acpc_WARPINV+orig.HEAD \
		../derivatives/"$sub"/"$img"-tbvu_avg_2acpc_warp.nii.gz \
		&> /dev/null

		rm -rf 	../derivatives/"$sub"/"$img"-tbvu_avg_2acpc_WARPINV+orig.*
		
		3dNwarpApply -nwarp  ../derivatives/"$sub"/"$img"-tbvu_avg_2acpc_warp.nii.gz \
		-source ../derivatives/"$sub"/"$img"-tbv_avg_2acpc.nii.gz \
		-prefix ../derivatives/"$sub"/"$img"-tbvu_avg_2acpc.nii.gz \
		-interp wsinc5 \
		&> /dev/null

		echo "$img: transform epi to template space"
		3dNwarpApply -nwarp  "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz \
		../derivatives/"$sub"/"$img"-tbvu_avg_2acpc_warp.nii.gz \
		../derivatives/"$sub"/"$img"-tbv_avg_2acpc.aff12.1D \
		../derivatives/"$sub"/"$img"-mocomatrix.aff12.1D" \
		-source ../derivatives/"$sub"/"$img"-tb.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-prefix ../derivatives/"$sub"/"$img"-tbvu_2mni.nii.gz \
		-dxyz "$epi_mni_res" \
		-interp wsinc5 \
		&> ../derivatives/"$sub"/"$img"-nonlin_epi2mni.log

		3dTstat \
		-prefix ../derivatives/"$sub"/"$img"-tbvu_2mni_avg.nii.gz \
		../derivatives/"$sub"/"$img"-tbvu_2mni.nii.gz \
		&> /dev/null

		echo "$img: transform epi mask to template space"
		3dNwarpApply -nwarp  "../derivatives/"$sub"/"$sub"_avg_T1w_2mni_brain_norm_WARP.nii.gz \
		../derivatives/"$sub"/"$img"-tbvu_avg_2acpc_warp.nii.gz \
		../derivatives/"$sub"/"$img"-tbv_avg_2acpc.aff12.1D" \
		-source ../derivatives/"$sub"/"$img"-tb_mask.nii.gz \
		-master "$ss_template_folder"/"$ss_template".nii.gz \
		-prefix ../derivatives/"$sub"/"$img"-tbvu_mask_2mni.nii.gz \
		-dxyz "$epi_mni_res" \
		-interp NN \
		&> ../derivatives/"$sub"/"$img"-nonlin_epimask2mni.log

		echo "$img: gaussian smoothing"
		3dBlurToFWHM \
		-input ../derivatives/"$sub"/"$img"-tbvu_2mni.nii.gz \
		-prefix ../derivatives/"$sub"/"$img"-tbvus_2mni.nii.gz \
		-mask ../derivatives/"$sub"/"$img"-tbvu_mask_2mni.nii.gz \
		-FWHM "$epi_fwhm" \
		-detrend \
		&> ../derivatives/"$sub"/"$img"-fwhm.log
		
		echo "$img: signal intensity scaling"
		3dTstat \
		-prefix ../derivatives/"$sub"/"$img"-tbvus_2mni_avg.nii.gz \
		../derivatives/"$sub"/"$img"-tbvus_2mni.nii.gz \
		&> /dev/null
		
		3dcalc -a ../derivatives/"$sub"/"$img"-tbvus_2mni.nii.gz \
		-b ../derivatives/"$sub"/"$img"-tbvus_2mni_avg.nii.gz \
		-expr '(a/b*100)-100' \
		-prefix ../derivatives/"$sub"/"$img"-tbvusn_2mni.nii.gz \
		&> ../derivatives/"$sub"/"$img"-signal_scaling.log
		
	done
	
	echo "$sub: create conjunction mask"
	3dmask_tool \
	-input ../derivatives/"$sub"/"$sub"*-tbvu_mask_2mni.nii.gz \
	-prefix ../derivatives/"$sub"/"$sub"_epi_brain_mask_2mni.nii.gz \
	-frac 1.0 \
	&> /dev/null

	echo "$sub: generate average brain extracted epi"
	3dMean -prefix ../derivatives/"$sub"/"$sub"_task-movie_allruns_bold-tbvu_2mni_avg.nii.gz \
	../derivatives/"$sub"/"$sub"_task-movie_run-??_bold-tbvu_2mni_avg.nii.gz \
	&> /dev/null

	cat ../derivatives/"$sub"/"$sub"*-mocopar.1D \
	> ../derivatives/"$sub"/"$sub"_epi-mocopar.1D

	echo "$sub: running deconvolve"
	3dDeconvolve \
	-jobs "$num_cpus" \
	-input ../derivatives/"$sub"/"$sub"*-tbvusn_2mni.nii.gz \
	-mask ../derivatives/"$sub"/"$sub"_epi_brain_mask_2mni.nii.gz \
	-polort A \
	-ortvec ../derivatives/"$sub"/"$sub"_epi-mocopar.1D no_interest \
	-errts ../derivatives/"$sub"/"$sub"_epi_cleaned_2mni.nii.gz \
	-nobucket \
	-x1D ../derivatives/"$sub"/"$sub"_epi_cleaned_2mni.xmat.1D \
	&> ../derivatives/"$sub"/"$sub"-deconvolve.log

	echo "$sub: running remlfit"	
	3dREMLfit -matrix ../derivatives/"$sub"/"$sub"_epi_cleaned_2mni.xmat.1D \
	-input "../derivatives/"$sub"/"$sub"*-tbvusn_2mni.nii.gz" \
	-mask ../derivatives/"$sub"/"$sub"_epi_brain_mask_2mni.nii.gz \
	-Rerrts ../derivatives/"$sub"/"$sub"_epi_cleaned_reml_2mni.nii.gz \
	&> ../derivatives/"$sub"/"$sub"-remlfit.log
	
done

echo "generate group average brain extracted anat"
3dMean -prefix ../derivatives/group_T1w_2mni_brain_norm.nii.gz \
../derivatives/sub-*/sub-*_avg_T1w_2mni_brain_norm.nii.gz \
&> /dev/null

echo "generate group average brain extracted epi"
3dMean -prefix ../derivatives/group_epi_2mni_brain.nii.gz \
../derivatives/sub-*/sub-*_task-movie_allruns_bold-tbvu_2mni_avg.nii.gz \
&> /dev/null

echo "generate group average anat"
3dMean -prefix ../derivatives/group_T1w_2mni_norm.nii.gz \
../derivatives/sub-*/sub-*_avg_T1w_2mni_norm.nii.gz \
&> /dev/null

echo "generate group average epi brain mask"
3dmask_tool -input ../derivatives/sub-*/sub-*_epi_brain_mask_2mni.nii.gz \
-prefix ../derivatives/group_epi_2mni_brain_mask.nii.gz \
-frac 1.0 \
&> /dev/null
