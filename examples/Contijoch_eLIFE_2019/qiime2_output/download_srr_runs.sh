# if there is any caching problem, run the following
# echo '/repository/user/main/public/root = "/tmp"' > $HOME/.ncbi/user-settings.mkfg
parallel --jobs 8 "fasterq-dump --split-files {}" ::: $(cat sra_ids.txt)
