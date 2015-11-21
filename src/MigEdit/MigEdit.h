#if !defined(__MigEdit_hdr__)
#define __MigEdit_hdr__

#include "MigEditIf.h"


class CMigEdit : public CMigEditIf
{
public:
    Error_t initInstance (/*enter parameters here*/);
    Error_t resetInstance ();

    //virtual Error_t process (float **ppfInputBuffer, float **ppfOutputBuffer, int iNumberOfFrames) = 0;

    CMigEdit ();
    virtual ~CMigEdit ();
};

#endif // #if !defined(__MigEdit_hdr__)



