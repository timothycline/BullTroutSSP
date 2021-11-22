#BullTrout SSP on Yeti

ssh tcline@yeti.cr.usgs.gov

#Create directory
mkdir ~/BullTroutSSP
#Write datafiles - Do in different terminal
scp -r Data tcline@yeti-dtn.cr.usgs.gov:~/BullTroutSSP


rm -r ~/BullTroutSSP
#clone git repository
cd ~/BullTroutSSP
ls
git clone https://github.com/timothycline/BullTroutSSP.git
git pull
git init -b main
git push --set-upstream origin main

