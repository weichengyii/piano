//#include "sndfile.h"
//
//PcmWriter::PcmWriter(char *filename, long size, int samplerate, int channels) {
//
//    total = 0;
//    info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
//    info.frames = size;
//    info.samplerate = samplerate;
//    info.channels = channels;
//    info.sections = 1;
//    info.seekable = 0;
//
//
//    if (!sf_format_check(&info))
//        info.format = (info.format & SF_FORMAT_TYPEMASK);
//
//    out = sf_open(filename, SFM_WRITE, &info);
//
//    if (!sf_format_check(&info))
//        perror("bad formatfor writing pcm");
//    if (!out) {
//        perror("cannot open file for writing");
//    }
//}
//
//long PcmWriter::write(double *data, long n) {
//    return sf_writef_double(out, (double *) data, n);
//}
//
//PcmWriter::~PcmWriter() {
//    close();
//}
//
//void PcmWriter::close() {
//    sf_close(out);
//}
//
//
