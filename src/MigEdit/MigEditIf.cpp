
// standard headers

// project headers
#include "MigEditConfig.h"

#include "ErrorDef.h"

#include "MigEditIf.h"
#include "MigEdit.h"

static const char*  kCMigEditBuildDate             = __DATE__;


//CMigEditIf::CMigEditIf ()
//{
//    // this never hurts
//    this->resetInstance ();
//}
//
//
//CMigEditIf::~CMigEditIf ()
//{
//    this->resetInstance ();
//}

const int  CMigEditIf::getVersion (const Version_t eVersionIdx)
{
    int iVersion = 0;

    switch (eVersionIdx)
    {
    case kMajor:
        iVersion    = MigEdit_VERSION_MAJOR; 
        break;
    case kMinor:
        iVersion    = MigEdit_VERSION_MINOR; 
        break;
    case kPatch:
        iVersion    = MigEdit_VERSION_PATCH; 
        break;
    case kNumVersionInts:
        iVersion    = -1;
        break;
    }

    return iVersion;
}
const char*  CMigEditIf::getBuildDate ()
{
    return kCMigEditBuildDate;
}

Error_t CMigEditIf::createInstance (CMigEditIf*& pCMigEdit)
{
    pCMigEdit = new CMigEdit ();

    if (!pCMigEdit)
        return kUnknownError;


    return kNoError;
}

Error_t CMigEditIf::destroyInstance (CMigEditIf*& pCMigEdit)
{
    if (!pCMigEdit)
        return kUnknownError;
    
    pCMigEdit->resetInstance ();
    
    delete pCMigEdit;
    pCMigEdit = 0;

    return kNoError;

}
