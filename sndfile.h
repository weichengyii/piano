#ifndef PCM_H
#define PCM_H

#include <sndfile.h>

class PcmWriter {
public:
    PcmWriter(char *filename, long, int, int);

    ~PcmWriter();

    long write(double *buf, long block_size);

    void close();

protected:
//    long total;
    SF_INFO info{};
    SNDFILE *out;

};

#endif
