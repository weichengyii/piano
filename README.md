Written by Clayton Otey, Dec 2007

Compilation Instructions:
First you need libsndfile.
Download the source from http://www.mega-nerd.com/libsndfile/ into a new directory.
cd <directory>
tar xzf *.tar.gz
./configure
make
sudo make install
This will install into /usr/local by default.  See the README and INSTALL files in the libsndfile source distribution for instructions on how to install on non-unix platforms or how to change the install locations.
If you're running linux and ld is not set up to look in /usr/local/lib, you may get an error like
"error while loading shared libraries libsndfile.so.1: cannot open shared object file: No such file or directory."
In this case you will need to add the following line to your /etc/ld.so.conf file (as root):
/usr/local/lib
and rebuild your ld.so cache by running 'sudo ldconfig'

Now, cd into the piano src directory and type `make'

Usage:
Run ./piano to see usage
When run, a mono .wav file called out.wav is created.  If the -b option is used a file called wave.out is generated instead, which is a binary file of doubleing point numbers in the range [-1:1].  This can be converted to any mono audio format using 3rd party tools e.g. libsndfile.