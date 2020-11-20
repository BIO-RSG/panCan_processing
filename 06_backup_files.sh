#! /bin/bash

# Copy files from one folder to another using the rsync command, which skips files that are identical between the two folders.

# Current options used: -auvz
# -a, --archive               archive mode; equals -rlptgoD (no -H,-A,-X)
    # -r, --recursive             recurse into directories
    # -l, --links                 copy symlinks as symlinks
    # -p, --perms                 preserve permissions
    # -t, --times                 preserve modification times
    # -g, --group                 preserve group
    # -o, --owner                 preserve owner (super-user only)
    # -D                          same as --devices --specials
        # --devices               preserve device files (super-user only)
        # --specials              preserve special files
# -u, --update                skip files that are newer on the receiver
# -v, --verbose               increase verbosity
# -z, --compress              compress file data during the transfer

in_path="/mnt/data3/claysa/"
out_path="/mnt/external/"


# Note (2020-11-12): If you are having issues with the backup, try changing the "auvz" option to "ruvz"


# A list of folders to backup, separated by spaces. Example formatting:
# to_backup=(climatologies/ OLCI-A/ OLCI-B/ MODIS/ SeaWiFS/ VIIRS-SNPP/)
to_backup=(climatologies/ OLCI-A/ OLCI-B/ MODIS/ SeaWiFS/ VIIRS-SNPP/)

for files in ${to_backup[@]}; do
	
	echo "Backing up $in_path$files to $out_path$files ..."
	rsync -ruvz $in_path$files $out_path$files
	
done
