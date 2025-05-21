Last login: Thu Feb 20 10:41:51 on ttys000
i439h@dkfz-vpn2067 Downloads % ssh i439h@odcf-worker02.dkfz.de        
i439h@odcf-worker02.dkfz.de's password: 
Last login: Wed Feb 19 11:01:08 2025 from b062-mb010-wl.inet.dkfz-heidelberg.de
####################################################################
#             ________ ______ __________________                   #
#             ___  __ \___  //_/___  ____/___  /                   #
#             __  / / /__  ,<   __  /_    __  /                    #
#             _  /_/ / _  /| |  _  __/    _  /__                   #
#             /_____/  /_/ |_|  /_/       /____/                   #
#       _______________                _____                       #
#       __  ____/___  /____  ____________  /______ ________        #
#       _  /     __  / _  / / /__  ___/_  __/_  _ \__  ___/        #
#       / /___   _  /  / /_/ / _(__  ) / /_  /  __/_  /            #
#       \____/   /_/   \__,_/  /____/  \__/  \___/ /_/             #
#                                                                  #
#                                                                  #
#                 Welcome to the ODCF Worker node                  #
#                                                                  #
#  ==============================================================  #
#   This host is a SHARED RESOURCE, mainly meant for development   #
#   interactive work and short term processing with low to mode-   #
#   rate resource requirements.                                    #
#                                                                  #
#   Long running processing and jobs requiring lots of CPU and     #
#   memory must be submitted AS REGULAR CLUSTER JOBS.              #
#                                                                  #
#   The cluster Wiki is at https://wiki.odcf.dkfz.de/pub/cluster/  #
#  ==============================================================  #
#                                                                  #
#  LOGIN TO THE CLUSTER IS ONLY GRANTED ON REQUEST.  Please        #
#  contact the email address below to request access. No           #
#  prerequisites are required.  Please mention your username.      #
#                                                                  #
#                                                                  #
#  If you have problems or need something installed, mail to:      #
#                cluster-support@dkfz-heidelberg.de                #
#                                                                  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#                   provided  and managed by the                   #
#        Omics IT and Data Management Core Facility (ODCF)         #
#                cluster-support@dkfz-heidelberg.de                #
####################################################################
-bash: warning: setlocale: LC_CTYPE: cannot change locale (UTF-8): No such file or directory


l
-bash-4.2$ 
-bash-4.2$ 
-bash-4.2$ l
total 64
drwxr-xr-x. 3 i439h B062  84 Aug  6  2023 nextflow
drwxr-xr-x. 2 i439h B062 382 Feb  6 11:25 projects
-bash-4.2$ l
total 64
drwxr-xr-x. 3 i439h B062  84 Aug  6  2023 nextflow
drwxr-xr-x. 2 i439h B062 382 Feb  6 11:25 projects
-bash-4.2$ cd projects/
-bash-4.2$ cd kitz_heroes/
bulkWGS/                      INFORM_HEROES_SNP_Genotyping/ singlecell_data/              
-bash-4.2$ cd kitz_heroes/singlecell_data/cellranger/
-bash-4.2$ l
total 2040
-rw-------. 1 i439h B062     86 Dec 12 14:36 650064[1].b06x-pbs01.ER
-rw-------. 1 i439h B062 161668 Dec 12 14:36 650064[1].b06x-pbs01.OU
-rw-------. 1 i439h B062    114 Dec 12 14:36 650064[2].b06x-pbs01.ER
-rw-------. 1 i439h B062 133588 Dec 12 14:35 650064[2].b06x-pbs01.OU
-rw-------. 1 i439h B062    114 Dec 12 14:36 650064[3].b06x-pbs01.ER
-rw-------. 1 i439h B062 114182 Dec 12 14:35 650064[3].b06x-pbs01.OU
-rw-------. 1 i439h B062      0 Dec 11 14:35 650064[4].b06x-pbs01.ER
-rw-------. 1 i439h B062  11501 Dec 11 14:35 650064[4].b06x-pbs01.OU
-rw-------. 1 i439h B062      0 Dec 11 14:35 650064[5].b06x-pbs01.ER
-rw-------. 1 i439h B062  11501 Dec 11 14:35 650064[5].b06x-pbs01.OU
-rw-------. 1 i439h B062    114 Jan 15 11:15 655176[1].b06x-pbs01.ER
-rw-------. 1 i439h B062 145875 Jan 15 11:15 655176[1].b06x-pbs01.OU
-rw-------. 1 i439h B062      0 Jan 14 11:14 655176[2].b06x-pbs01.ER
-rw-------. 1 i439h B062     39 Jan 14 11:14 655176[2].b06x-pbs01.OU
-rw-------. 1 i439h B062    114 Jan 22 11:17 658453[1].b06x-pbs01.ER
-rw-------. 1 i439h B062 113720 Jan 22 11:17 658453[1].b06x-pbs01.OU
-rw-------. 1 i439h B062      0 Jan 21 11:16 658453[2].b06x-pbs01.ER
-rw-------. 1 i439h B062 127883 Jan 21 15:45 658453[2].b06x-pbs01.OU
-rw-------. 1 i439h B062   1784 Dec 11 14:35 __pool11.mro
-rw-------. 1 i439h B062   1784 Dec 11 14:35 __pool19.mro
-rw-------. 1 i439h B062   1784 Dec 11 14:35 __pool23.mro
-rw-------. 1 i439h B062   1784 Jan 21 11:16 __pool27.mro
-rw-------. 1 i439h B062      0 Nov 18 11:47 cellranger_error_%J_%I.err
-rw-------. 1 i439h B062   7377 Nov 18 11:49 cellranger_output_%J_%I.out
-rwxr-xr-x. 1 i439h B062     56 Nov 18 11:56 job.sh
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool10.csv
drwx--S---. 5 i439h B062    292 Dec 12 14:36 pool11
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool11.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool12.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool13.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool14.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool15.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool16.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool17.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool18.csv
drwx--S---. 5 i439h B062    292 Dec 12 14:35 pool19
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool19.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool20.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool21.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool22.csv
drwx--S---. 5 i439h B062    292 Dec 12 14:35 pool23
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool23.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool24.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool25.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool26.csv
drwx--S---. 5 i439h B062    292 Jan 22 11:17 pool27
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool27.csv
drwx--S---. 4 i439h B062    532 Dec  9 18:34 pool28
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool28.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool29.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool30.csv
drwx--S---. 4 i439h B062    494 Dec  9 10:33 pool31
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool31.csv
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool32.csv
drwx--S---. 4 i439h B062    518 Jan 16 16:22 pool33
-rwxr-xr-x. 1 i439h B062    314 Nov 18 11:07 pool33.csv
-rwxr-xr-x. 1 i439h B062    286 Nov 18 11:07 pool34.csv
-rwxr-xr-x. 1 i439h B062    286 Nov 18 11:07 pool35.csv
drwx--S---. 4 i439h B062    532 Dec 11 14:36 pool36
-rwxr-xr-x. 1 i439h B062    286 Nov 18 11:07 pool36.csv
drwx--S---. 4 i439h B062    494 Dec 11 14:36 pool37
-rwxr-xr-x. 1 i439h B062    286 Nov 18 11:07 pool37.csv
drwx--S---. 4 i439h B062    532 Nov 29 23:38 pool38
-rwxr-xr-x. 1 i439h B062    286 Nov 18 11:07 pool38.csv
drwx--S---. 4 i439h B062    532 Nov 30 00:21 pool39
-rwxr-xr-x. 1 i439h B062    286 Nov 18 11:07 pool39.csv
-rwxr-xr-x. 1 i439h B062    286 Nov 18 11:07 pool40.csv
drwx--S---. 4 i439h B062    532 Nov 30 14:31 pool41
-rwxr-xr-x. 1 i439h B062    286 Nov 18 11:07 pool41.csv
-rwxr-xr-x. 1 i439h B062    299 Nov 18 11:07 pool5.csv
-rwxr-xr-x. 1 i439h B062    298 Nov 18 11:07 pool6.csv
-rwxr-xr-x. 1 i439h B062    298 Nov 18 11:07 pool7.csv
-rwxr-xr-x. 1 i439h B062    298 Nov 18 11:07 pool8.csv
-rwxr-xr-x. 1 i439h B062    310 Nov 18 11:07 pool9.csv
-rwxr--r--. 1 i439h B062   1472 Jan 21 11:15 run_cellranger_arc_kitz_cluster.sh
-bash-4.2$ 
-bash-4.2$ cd ..
-bash-4.2$ l
total 2344
-rw-rw-r--.  1 s472a B062  15984 Sep 11  2021 6909_JW_WES_HighPass_exome.xlsx
drwxr-sr-x.  3 s472a B062     30 Nov 12 13:40 AG_Thongjuea
drwxrwsr-x. 13 s472a B062    597 Aug 28  2023 BioSkryb
-rw-rw-r--.  1 s472a B062   5985 Sep 11  2021 FASTQtransfer_sarcoma_Pfister_20210911_lowpass_forCNV.csv
drwxrwsr-x. 14 i439h B062   2342 Jan 22 16:19 cellranger
-rw-rw-r--.  1 s472a B062 775553 Sep 11  2021 sarcoma_20210911.pdf
-rw-rw-r--.  1 s472a B062 793055 Sep 13  2021 sarcoma_pilot_study_data_transfer_and_initial_analysis.zip
drwxr-sr-x.  6 s472a B062    220 Nov  1  2023 scWES
drwxrwsr-x.  2 s472a B062   1161 Nov 13 10:43 tools
-bash-4.2$ chmod 775 -R cellranger
-bash-4.2$ client_loop: send disconnect: Broken pipe
i439h@dkfz-vpn2067 Downloads % 
  [Restored 9. Mar 2025 at 16:31:43]
Last login: Sun Mar  9 16:31:42 on ttys002
i439h@b062-mb010-wl Downloads % 
  [Restored 7. Apr 2025 at 14:37:25]
Last login: Mon Apr  7 14:37:25 on ttys000
i439h@b062-mb010-wl Downloads % 
  [Restored 11. Apr 2025 at 22:51:57]
Last login: Fri Apr 11 22:51:57 on ttys000
i439h@b062-mb010-wl Downloads % 
  [Restored 2. May 2025 at 12:55:22]
Last login: Fri May  2 12:55:22 on ttys000
i439h@b062-mb010-wl Downloads % 
  [Restored 18. May 2025 at 17:46:05]
Last login: Sun May 18 17:46:02 on ttys000
i439h@b062-mb010-wl Downloads % ssh i439h@odcf-worker02.dkfz.de        
i439h@odcf-worker02.dkfz.de's password: 
Last login: Wed May 21 12:22:35 2025 from dkfz-vpn2045.inet.dkfz-heidelberg.de
####################################################################
#             ________ ______ __________________                   #
#             ___  __ \___  //_/___  ____/___  /                   #
#             __  / / /__  ,<   __  /_    __  /                    #
#             _  /_/ / _  /| |  _  __/    _  /__                   #
#             /_____/  /_/ |_|  /_/       /____/                   #
#       _______________                _____                       #
#       __  ____/___  /____  ____________  /______ ________        #
#       _  /     __  / _  / / /__  ___/_  __/_  _ \__  ___/        #
#       / /___   _  /  / /_/ / _(__  ) / /_  /  __/_  /            #
#       \____/   /_/   \__,_/  /____/  \__/  \___/ /_/             #
#                                                                  #
#                                                                  #
#                 Welcome to the ODCF Worker node                  #
#                                                                  #
#  ==============================================================  #
#   This host is a SHARED RESOURCE, mainly meant for development   #
#   interactive work and short term processing with low to mode-   #
#   rate resource requirements.                                    #
#                                                                  #
#   Long running processing and jobs requiring lots of CPU and     #
#   memory must be submitted AS REGULAR CLUSTER JOBS.              #
#                                                                  #
#   The cluster Wiki is at https://wiki.odcf.dkfz.de/pub/cluster/  #
#  ==============================================================  #
#                                                                  #
#  LOGIN TO THE CLUSTER IS ONLY GRANTED ON REQUEST.  Please        #
#  contact the email address below to request access. No           #
#  prerequisites are required.  Please mention your username.      #
#                                                                  #
#                                                                  #
#  If you have problems or need something installed, mail to:      #
#                cluster-support@dkfz-heidelberg.de                #
#                                                                  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
#                   provided  and managed by the                   #
#        Omics IT and Data Management Core Facility (ODCF)         #
#                cluster-support@dkfz-heidelberg.de                #
####################################################################
-bash: warning: setlocale: LC_CTYPE: cannot change locale (UTF-8): No such file or directory
-bash-4.2$ 
-bash-4.2$ 
-bash-4.2$ cd projects/Emseq_temp/
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED  67 May 21 12:29 README.md
drwxr-sr-x. 12 i439h W610-B062-PED 341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED  27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED 137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED 214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED 230 Jan 23 12:36 tools
-bash-4.2$ l codes
total 424
drwxr-sr-x. 2 i439h W610-B062-PED  178 Nov 19  2024 biscuit_pipeline
drwxr-sr-x. 7 i439h W610-B062-PED  122 Nov 15  2024 cnv_calling
-rwxr-xr-x. 1 i439h W610-B062-PED  205 Oct 10  2024 codes.Rproj
drwxr-sr-x. 2 i439h W610-B062-PED  463 Jan 27 13:37 dmr
-rwxr-xr-x. 1 i439h W610-B062-PED 7167 Oct 10  2024 dmr.R
drwxr-sr-x. 2 i439h W610-B062-PED   41 Feb  4 16:38 integration
drwxr-sr-x. 2 i439h W610-B062-PED   32 Oct 10  2024 methylation_array
drwxr-sr-x. 4 i439h W610-B062-PED  432 Dec  2 20:07 nextflow
drwxr-sr-x. 2 i439h W610-B062-PED   33 Nov 29 11:08 qc_summary
drwxr-sr-x. 2 i439h W610-B062-PED   75 Oct 23  2024 sample_prep
drwxr-sr-x. 3 i439h W610-B062-PED   25 Nov 15  2024 snv_calling
-bash-4.2$ 
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED  67 May 21 12:29 README.md
drwxr-sr-x. 12 i439h W610-B062-PED 341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED  27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED 137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED 214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED 230 Jan 23 12:36 tools
-bash-4.2$ pwd
/home/i439h/projects/Emseq_temp
-bash-4.2$ cd ..
-bash-4.2$ l
total 336
lrwxrwxrwx. 1 i439h B062 88 Sep 20  2023 EMseq -> /omics/odcf/analysis/OE0290_projects/pediatric_tumor/whole_genome_biosulfite_sequencing/
lrwxrwxrwx. 1 i439h B062 64 Oct 15  2024 Emseq_temp -> /omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/
lrwxrwxrwx. 1 i439h B062 30 Nov 18  2024 HEROES-AYA -> /b06x-isi/b062/g-i/HEROES-AYA/
lrwxrwxrwx. 1 i439h B062 25 Nov 27 15:12 cellMeth -> Emseq_temp/tools/cellMeth
lrwxrwxrwx. 1 i439h B062 88 May 29  2024 demultiplex_pipeline -> /home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline
lrwxrwxrwx. 1 i439h B062 44 Nov 22 11:36 demuxit -> /home/i439h/projects/pool_temp/tools/demuxit
lrwxrwxrwx. 1 i439h B062 47 Sep 18  2023 heroes-aya -> /omics/odcf/analysis/OE0290_projects/heroes-aya
lrwxrwxrwx. 1 i439h B062 78 Dec  2 14:09 heroes_temp -> /omics/odcf/analysis/OE0290_projects_temp/heroes-aya/AG_Thongjuea/Sub_project1
lrwxrwxrwx. 1 i439h B062 37 May 13  2024 hipo_k35 -> /omics/odcf/analysis/hipo2/hipo_K35R/
lrwxrwxrwx. 1 i439h B062 42 Oct 15  2024 hipo_temp -> /omics/odcf/analysis/hipo2_temp/hipo_K35R/
lrwxrwxrwx. 1 i439h B062 30 Nov 18  2024 kitz_heroes -> /b06x-isi/b062/g-i/HEROES-AYA/
lrwxrwxrwx. 1 i439h B062 67 Oct 15  2024 pool_temp -> /omics/odcf/analysis/OE0290_projects_temp/heroes-aya_pools/sadeghi/
lrwxrwxrwx. 1 i439h B062 54 Sep 18  2023 pools -> /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/
lrwxrwxrwx. 1 i439h B062 49 Apr 17  2024 sarcoma -> /omics/odcf/analysis/OE0290_projects/sarcoma/Iman
-bash-4.2$ l 
total 336
lrwxrwxrwx. 1 i439h B062 88 Sep 20  2023 EMseq -> /omics/odcf/analysis/OE0290_projects/pediatric_tumor/whole_genome_biosulfite_sequencing/
lrwxrwxrwx. 1 i439h B062 64 Oct 15  2024 Emseq_temp -> /omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/
lrwxrwxrwx. 1 i439h B062 30 Nov 18  2024 HEROES-AYA -> /b06x-isi/b062/g-i/HEROES-AYA/
lrwxrwxrwx. 1 i439h B062 25 Nov 27 15:12 cellMeth -> Emseq_temp/tools/cellMeth
lrwxrwxrwx. 1 i439h B062 88 May 29  2024 demultiplex_pipeline -> /home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline
lrwxrwxrwx. 1 i439h B062 44 Nov 22 11:36 demuxit -> /home/i439h/projects/pool_temp/tools/demuxit
lrwxrwxrwx. 1 i439h B062 47 Sep 18  2023 heroes-aya -> /omics/odcf/analysis/OE0290_projects/heroes-aya
lrwxrwxrwx. 1 i439h B062 78 Dec  2 14:09 heroes_temp -> /omics/odcf/analysis/OE0290_projects_temp/heroes-aya/AG_Thongjuea/Sub_project1
lrwxrwxrwx. 1 i439h B062 37 May 13  2024 hipo_k35 -> /omics/odcf/analysis/hipo2/hipo_K35R/
lrwxrwxrwx. 1 i439h B062 42 Oct 15  2024 hipo_temp -> /omics/odcf/analysis/hipo2_temp/hipo_K35R/
lrwxrwxrwx. 1 i439h B062 30 Nov 18  2024 kitz_heroes -> /b06x-isi/b062/g-i/HEROES-AYA/
lrwxrwxrwx. 1 i439h B062 67 Oct 15  2024 pool_temp -> /omics/odcf/analysis/OE0290_projects_temp/heroes-aya_pools/sadeghi/
lrwxrwxrwx. 1 i439h B062 54 Sep 18  2023 pools -> /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/
lrwxrwxrwx. 1 i439h B062 49 Apr 17  2024 sarcoma -> /omics/odcf/analysis/OE0290_projects/sarcoma/Iman
-bash-4.2$ cd Emseq_temp/
-bash-4.2$ k
-bash: k: command not found
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED  67 May 21 12:29 README.md
drwxr-sr-x. 12 i439h W610-B062-PED 341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED  27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED 137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED 214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED 230 Jan 23 12:36 tools
-bash-4.2$ l codes/
biscuit_pipeline/  codes.Rproj        dmr.R              methylation_array/ qc_summary/        sample_prep/       
cnv_calling/       dmr/               integration/       nextflow/          .Rproj.user/       snv_calling/       
-bash-4.2$ l codes/cnv_calling/
total 200
drwxr-sr-x. 2 i439h W610-B062-PED 31 Oct 10  2024 QDNASeq
drwxr-sr-x. 2 i439h W610-B062-PED 54 Jan  8 21:53 cfdna
drwxr-sr-x. 2 i439h W610-B062-PED 30 Oct 10  2024 cnMops
drwxr-sr-x. 2 i439h W610-B062-PED 52 Nov 15  2024 gistic2
drwxr-sr-x. 2 i439h W610-B062-PED 68 Oct 10  2024 hatchet
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED  67 May 21 12:29 README.md
drwxr-sr-x. 12 i439h W610-B062-PED 341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED  27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED 137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED 214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED 230 Jan 23 12:36 tools
-bash-4.2$ l data
total 40
drwxr-sr-x. 2 i439h W610-B062-PED 48 Jan 27 13:35 processed
-bash-4.2$ l data/processed/
total 140928
-rw-r--r--. 1 i439h W610-B062-PED 116616971 Jan 28 13:27 normalized_methyl_region.rdata
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED  67 May 21 12:29 README.md
drwxr-sr-x. 12 i439h W610-B062-PED 341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED  27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED 137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED 214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED 230 Jan 23 12:36 tools
-bash-4.2$ l datasets/
total 152
-rwxr-xr-x. 1 i439h W610-B062-PED 8351 Oct 23  2024 HEROES_AYA_LIQUID_BIOPSY.csv
drwxr-sr-x. 2 i439h W610-B062-PED   39 Nov 27 15:32 masked_regions
drwxr-sr-x. 2 i439h W610-B062-PED    0 Nov 29 10:58 references
drwxr-sr-x. 2 i439h W610-B062-PED   74 Dec  2 11:41 sample_sheets
-bash-4.2$ l datasets/sample_sheets/
total 64
-rwxrwxr-x. 1 i439h W610-B062-PED 4529 Dec  2 11:54 sampleSheet_undone.csv
-rwxrwxr-x. 1 i439h W610-B062-PED 2784 Nov 11  2024 sample_sheet.csv
-bash-4.2$ l datasets/references/
total 0
-bash-4.2$ l datasets/masked_regions/
total 64
-rw-r--r--. 1 i439h W610-B062-PED 26898 Nov 27 15:31 hg38-blacklist.v2.bed
-bash-4.2$ l 
total 232
-rw-r--r--.  1 i439h W610-B062-PED  67 May 21 12:29 README.md
drwxr-sr-x. 12 i439h W610-B062-PED 341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED  27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED 137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED 214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED 230 Jan 23 12:36 tools
-bash-4.2$ l results/
total 320
drwxr-sr-x. 3 i439h W610-B062-PED  29 Jan 23 11:52 QC
drwxr-sr-x. 5 i439h W610-B062-PED  76 Nov 23 00:44 biscuit_pipeline
drwxr-sr-x. 7 i439h W610-B062-PED 128 Jan  2 13:04 cnv_calling
drwxr-sr-x. 2 i439h W610-B062-PED  34 Jan 27 10:53 dmr
drwxr-sr-x. 3 i439h W610-B062-PED  21 Jan 27 13:55 figures
drwxr-sr-x. 3 i439h W610-B062-PED  37 Feb  4 16:24 methylation_array
drwxr-sr-x. 8 i439h W610-B062-PED 159 Oct 23  2024 nextflow
drwxr-sr-x. 3 i439h W610-B062-PED  21 Jan 27 13:56 tables
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED  67 May 21 12:29 README.md
drwxr-sr-x. 12 i439h W610-B062-PED 341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED  27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED 137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED 214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED 230 Jan 23 12:36 tools
-bash-4.2$ l results/nextflow/
total 240
drwxr-sr-x.  8 i439h W610-B062-PED  167 Oct 23  2024 bismark
drwxr-sr-x.  3 i439h W610-B062-PED 2214 Dec  3 02:08 fastqc
drwxr-sr-x.  3 i439h W610-B062-PED   25 Oct 23  2024 multiqc
drwxr-sr-x.  2 i439h W610-B062-PED 1930 Dec 13 00:13 pipeline_info
drwxr-sr-x. 21 i439h W610-B062-PED  830 Dec 11 10:00 qualimap
drwxr-sr-x.  4 i439h W610-B062-PED   46 Oct 15  2024 trimgalore
-bash-4.2$ l results/nextflow/bismark/
total 240
drwxr-sr-x. 3 i439h W610-B062-PED   22 Oct 16  2024 alignments
drwxr-sr-x. 3 i439h W610-B062-PED 1612 Dec 11 05:01 deduplicated
drwxr-sr-x. 7 i439h W610-B062-PED  156 Oct 22  2024 methylation_calls
drwxr-sr-x. 3 i439h W610-B062-PED  934 Dec 12 03:48 preseq
drwxr-sr-x. 2 i439h W610-B062-PED 1495 Dec 13 00:03 reports
drwxr-sr-x. 2 i439h W610-B062-PED   89 Dec 13 00:03 summary
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED  67 May 21 12:29 README.md
drwxr-sr-x. 12 i439h W610-B062-PED 341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED  27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED 137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED 214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED 230 Jan 23 12:36 tools
-bash-4.2$ l results/
total 320
drwxr-sr-x. 3 i439h W610-B062-PED  29 Jan 23 11:52 QC
drwxr-sr-x. 5 i439h W610-B062-PED  76 Nov 23 00:44 biscuit_pipeline
drwxr-sr-x. 7 i439h W610-B062-PED 128 Jan  2 13:04 cnv_calling
drwxr-sr-x. 2 i439h W610-B062-PED  34 Jan 27 10:53 dmr
drwxr-sr-x. 3 i439h W610-B062-PED  21 Jan 27 13:55 figures
drwxr-sr-x. 3 i439h W610-B062-PED  37 Feb  4 16:24 methylation_array
drwxr-sr-x. 8 i439h W610-B062-PED 159 Oct 23  2024 nextflow
drwxr-sr-x. 3 i439h W610-B062-PED  21 Jan 27 13:56 tables
-bash-4.2$ l results/figures/
total 40
drwxr-sr-x. 2 i439h W610-B062-PED 250 Feb  3 14:28 dmr
-bash-4.2$ l results/figures/dmr/
total 8648
-rw-r--r--. 1 i439h W610-B062-PED    5350 Jan 28 15:12 annotation_dmr_reg_plot.pdf
-rw-r--r--. 1 i439h W610-B062-PED    5250 Jan 28 14:52 diffMethPerChr_region.pdf
-rw-r--r--. 1 i439h W610-B062-PED 6865347 Jan 28 16:28 dmr_reg_heatmap.pdf
-rw-r--r--. 1 i439h W610-B062-PED    5009 Jan 28 15:11 dmrs_annotation_summary.pdf
-rw-r--r--. 1 i439h W610-B062-PED     491 Jan 28 15:13 dmrs_reg_detailed_annotation.csv
-rw-r--r--. 1 i439h W610-B062-PED    8226 Feb  3 14:28 pca_plot.pdf
-bash-4.2$ l results/tables/
total 40
drwxr-sr-x. 2 i439h W610-B062-PED 108 Jan 27 13:57 dmr
-bash-4.2$ l results/tables/dmr/
total 8360
-rw-r--r--. 1 i439h W610-B062-PED 2921723 Jan 28 14:52 DMRs_region.txt
-rw-r--r--. 1 i439h W610-B062-PED 1130373 Jan 28 15:11 dmr_reg_obj.Rds
-rw-r--r--. 1 i439h W610-B062-PED 3667232 Jan 28 15:11 dmr_region_annotated.Rds
-bash-4.2$ pwd
/home/i439h/projects/Emseq_temp
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED  67 May 21 12:29 README.md
drwxr-sr-x. 12 i439h W610-B062-PED 341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED  27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED 137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED 214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED 230 Jan 23 12:36 tools
-bash-4.2$ nano README.md 
-bash-4.2$ 
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED 1536 May 21 14:53 README.md
drwxr-sr-x. 12 i439h W610-B062-PED  341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED   27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED  137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED  214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED  230 Jan 23 12:36 tools
-bash-4.2$ less README.md 
-bash-4.2$ cat > README.md << 'EOF'
> # HEROES_EMseq
> 
> **EMSeq Data Analysis for Liquid Biopsy in Pediatric Tumors**
> 
> This repository contains code and results for analysis of EMSeq (Enzymatic Methyl-seq) data from liquid biopsy samples in pediatric tumors.
> 
> ---
> 
> ## Table of Contents
> 1. [Project Overview](#project-overview)  
> 2. [Directory Structure](#directory-structure)  
> 3. [Code Modules](#code-modules)  
> 4. [Data Files](#data-files)  
> 5. [Results](#results)  
> 6. [Usage](#usage)  
> 7. [Contributing](#contributing)  
> 8. [License](#license)  
> 
> ---
> 
> ## Project Overview
> The goal of this project is to perform comprehensive DNA methylation analysis using EMSeq data derived from liquid biopsy samples of pediatric tumors. Analyses include read alignment, quality control, differential methylation region (DMR) detection, and copy number variation (CNV) calling.
> 
> **Main directories:**
> - `/omics/odcf/analysis/OE0290_projects/pediatric_tumor/whole_genome_biosulfite_sequencing/`
> - Temporary working directory with up-to-date code and results: `/omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/`
> 
> ---
> 
> ## Directory Structure
> \`\`\`
> ├── codes/                  # Analysis pipelines and scripts
> │   ├── nextflow/           # nf-core/methylseq pipeline
> │   ├── dmr/                # Differentially Methylated Region analysis
> │   ├── qc_summary/         # Sample-level QC summary scripts
> │   ├── sample_prep/        # Input preparation for pipelines
> │   └── cnv_calling/        # CNV calling with cfDNA-specific tools
> │       └── cfdna/          # Main CNV-calling tool (preferred for liquid biopsy)
> │
> ├── data/                   # Preprocessed methylation calls (normalized/loading)
> │   └── methylation_calls/  # Input for DMR analysis
> │
> ├── datasets/               # Sample sheets for pipeline inputs
> │   └── samples_sheet.tsv   # Metadata and sample definitions
> │
> ├── results/                # Analysis outputs
> │   ├── nextflow/           # nf-core pipeline outputs
> │   │   ├── bismark/        # Alignment and methylation calls
> │   │   │   ├── deduplicated/   # BAM files for CNV calling
> │   │   │   └── methylation_calls/  # Files for DMR analysis
> │   │   └── ...
> │   ├── dmr/                # DMR analysis outputs and figures
> │   │   └── figures/        # Visualization of DMR results
> │   └── biscuit_pipeline/   # Test scripts for Biscuit-based pipeline
> │
> └── README.md               # This file
> \`\`\`
> 
> ---
> 
> ## Code Modules
> 
> ### 1. \`codes/nextflow\`
> - Implements the [nf-core/methylseq](https://github.com/nf-core/methylseq) pipeline using Nextflow.
> - Produces alignment (Bismark), deduplication, and methylation call outputs.
> 
> ### 2. \`codes/dmr\`
> - Scripts for identifying and annotating differentially methylated regions across samples.
> 
> ### 3. \`codes/qc_summary\`
> - Aggregates QC metrics (e.g., coverage, duplication rate) across EMSeq samples into summary reports.
> 
> ### 4. \`codes/sample_prep\`
> - Generates input manifests and FASTQ file lists required by the Nextflow pipeline.
> 
> ### 5. \`codes/cnv_calling\`
> - Focuses on CNV detection from cfDNA data.
> - **Use the \`cfdna/\` tool**—other scripts or tools in this directory are deprecated for liquid biopsy analyses.
> 
> ---
> 
> ## Data Files
> - **Normalized methylation calls** are stored in \`data/methylation_calls/\`. These tables are used as input for downstream DMR analysis in the \`codes/dmr/\` module.
> - **Sample sheets** in \`datasets/\` define sample metadata (IDs, paths, groups) for pipeline runs.
> 
> ---
> 
> ## Results
> 
> ### Nextflow Pipeline Outputs (\`results/nextflow\`)
> - **\`bismark/deduplicated/\`**: Deduplicated BAM files for CNV calling.
> - **\`bismark/methylation_calls/\`**: Cytosine call files for DMR analysis.
> 
> ### DMR Analysis (\`results/dmr\`)
> - **Figures**: Plots and heatmaps summarizing differentially methylated regions.
> 
> ### Biscuit Pipeline (\`results/biscuit_pipeline\`)
> - Prototype scripts for running the [Biscuit](https://informatics.fas.harvard.edu/biscuit/) toolchain. Modify and re-run these to benchmark against Nextflow results.
> 
> ---
> 
> ## Usage
> 1. **Prepare samples**: Edit \`datasets/samples_sheet.tsv\` with sample metadata.  
> 2. **Run Nextflow**:  
>    \`\`\`bash
>    cd codes/nextflow
>    nextflow run nf-core/methylseq -profile odcf --input ../../datasets/samples_sheet.tsv
>    \`\`\`
> 3. **QC Summary**:  
>    \`\`\`bash
>    cd codes/qc_summary
>    Rscript summarize_qc.R --input ../../results/nextflow/bismark/
>    \`\`\`
> 4. **DMR Analysis**:  
>    \`\`\`bash
>    cd codes/dmr
>    python run_dmr.py --calls ../../data/methylation_calls/ --output ../../results/dmr/
>    \`\`\`
> 5. **CNV Calling**:  
>    \`\`\`bash
>    cd codes/cnv_calling/cfdna
>    bash run_cfdna.sh --bam ../../results/nextflow/bismark/deduplicated/*.bam
>    \`\`\`
> 
> ---
> 
> ## Contributing
> - Please open issues or submit pull requests for improvements or bug fixes.
> - Follow [ODCF coding standards](https://odcf.org/coding-standards).
> 
> ---
> 
> ## License
> This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
> EOF
-bash-4.2$ l
total 232
-rw-r--r--.  1 i439h W610-B062-PED 5150 May 21 15:12 README.md
drwxr-sr-x. 12 i439h W610-B062-PED  341 Feb  4 16:37 codes
drwxr-sr-x.  3 i439h W610-B062-PED   27 Jan 27 13:35 data
drwxr-sr-x.  5 i439h W610-B062-PED  137 Nov 29 10:50 datasets
drwxr-sr-x. 10 i439h W610-B062-PED  214 Feb  4 16:03 results
drwxr-sr-x. 10 i439h W610-B062-PED  230 Jan 23 12:36 tools
-bash-4.2$ less README.md 
-bash-4.2$ less README.md 

- Aggregates QC metrics (e.g., coverage, duplication rate) across EMSeq samples into summary reports.

### 4. \`codes/sample_prep\`
- Generates input manifests and FASTQ file lists required by the Nextflow pipeline.

### 5. \`codes/cnv_calling\`
- Focuses on CNV detection from cfDNA data.
- **Use the \`cfdna/\` tool**<E2><80><94>other scripts or tools in this directory are deprecated for liquid biopsy analyses.

---

## Data Files
- **Normalized methylation calls** are stored in \`data/methylation_calls/\`. These tables are used as input for downstream DMR analysis in the \`codes/dmr/\` module.
- **Sample sheets** in \`datasets/\` define sample metadata (IDs, paths, groups) for pipeline runs.

---

## Results

### Nextflow Pipeline Outputs (\`results/nextflow\`)
- **\`bismark/deduplicated/\`**: Deduplicated BAM files for CNV calling.
- **\`bismark/methylation_calls/\`**: Cytosine call files for DMR analysis.

### DMR Analysis (\`results/dmr\`)
- **Figures**: Plots and heatmaps summarizing differentially methylated regions.

### Biscuit Pipeline (\`results/biscuit_pipeline\`)
- Prototype scripts for running the [Biscuit](https://informatics.fas.harvard.edu/biscuit/) toolchain. Modify and re-run these to benchmark against Nextflow results.

---

## Usage
1. **Prepare samples**: Edit \`datasets/samples_sheet.tsv\` with sample metadata.  
2. **Run Nextflow**:  
   \`\`\`bash
   cd codes/nextflow
   nextflow run nf-core/methylseq -profile odcf --input ../../datasets/samples_sheet.tsv
   \`\`\`
3. **QC Summary**:  
   \`\`\`bash
   cd codes/qc_summary
   Rscript summarize_qc.R --input ../../results/nextflow/bismark/
   \`\`\`
4. **DMR Analysis**:  
   \`\`\`bash
   cd codes/dmr
   python run_dmr.py --calls ../../data/methylation_calls/ --output ../../results/dmr/
   \`\`\`
5. **CNV Calling**:  
   \`\`\`bash
   cd codes/cnv_calling/cfdna
   bash run_cfdna.sh --bam ../../results/nextflow/bismark/deduplicated/*.bam
   \`\`\`

---

## Contributing
- Please open issues or submit pull requests for improvements or bug fixes.
- Follow [ODCF coding standards](https://odcf.org/coding-standards).

---

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
