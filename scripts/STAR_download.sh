#!/bin/bash

# Move to home directory
cd /home/bhavana || exit

# Download STAR 2.7.10b
wget -P /home/bhavana https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10b.tar.gz

# Extract the archive
tar -xvzf 2.7.10b.tar.gz

# Compile STAR
cd STAR-2.7.10b/source || exit
make STAR

# Add STAR to PATH (if not already present)
if ! grep -q 'STAR-2.7.10b/source' ~/.bashrc; then
    echo 'export PATH=/home/bhavana/STAR-2.7.10b/source:$PATH' >> ~/.bashrc
    echo "✅ Added STAR to your PATH in ~/.bashrc"
else
    echo "ℹ️  STAR already in your PATH"
fi

# Reload shell settings
source ~/.bashrc

# Confirmation message
echo "✅ STAR has been successfully downloaded, compiled, and added to your PATH."
