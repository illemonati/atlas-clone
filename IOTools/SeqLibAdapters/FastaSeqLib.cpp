#include "FastaSeqLib.h"


bool FastaSeqLib::GetSequence(const int &refid, const int &start, const int &stop, std::string &ref)
{
    //reference.Rewind();
    if(last_buffer_seq>=refid){
        //look in the buffer
        ref= buffer.at(refid).substr(start, (stop-start+1));
        return true;
    }else{
        //fetch new sequence
        SeqLib::UnalignedSequence seq;
        bool res;
        do{
            res=reference.GetNextSequence(seq);
            ref=seq.Seq.substr(start, (stop-start+1));
            buffer.push_back(seq.Seq);
            last_buffer_seq++;
        }while(refid<last_buffer_seq && res);
        return res;
    }
}

bool FastaSeqLib::Open(std::string &filepath, std::string &index)
{
    return reference.Open(filepath);
}

bool FastaSeqLib::GetBase(const int &refid, const int &position, char &base)
{

    if(last_buffer_seq>=refid){
        //look in the buffer
        base= buffer.at(refid).at(position);
        return true;
    }else{
        //fetch new sequence
        SeqLib::UnalignedSequence seq;
        bool res;
        do{
            res=reference.GetNextSequence(seq);
            base=seq.Seq.at(position);
            buffer.push_back(seq.Seq);
            last_buffer_seq++;
        }while(refid<last_buffer_seq && res);
        return res;
    }
}
