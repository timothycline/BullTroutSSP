#BullTrout SSP on Yeti

ssh tcline@yeti.cr.usgs.gov

#Create directory
mkdir ~/BullTroutSSP
#Write datafiles - Do in different terminal
scp -r Data tcline@yeti-dtn.cr.usgs.gov:~/BullTroutSSP

#clone git repository
cd ~/BullTroutSSP

git clone https://github.com/timothycline/BullTroutSSP.git
