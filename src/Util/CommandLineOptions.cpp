
#include "Util.h"
#include "CommandLineOptions.h"

CCommandLineOptions::CCommandLineOptions() :
    m_iNumOfOptions(0),
    m_ppCOptions(0),
    m_ppCClArgs(0)
{
    resetInstance ();
}

Error_t CCommandLineOptions::createInstance( CCommandLineOptions*& pCCommandLineOptions )
{
    pCCommandLineOptions = new CCommandLineOptions ();

    if (!pCCommandLineOptions)
        return kMemError;

    return kNoError;
}

Error_t CCommandLineOptions::destroyInstance( CCommandLineOptions*& pCCommandLineOptions )
{
    if (!pCCommandLineOptions)
        return kNoError;

    delete pCCommandLineOptions;
    pCCommandLineOptions   = 0;

    return kNoError;
}

Error_t CCommandLineOptions::initInstance (const COption myOptions[], int iNumOfOptions)
{
    Error_t  rErr = kNoError;

    this->resetInstance();

    m_ppCOptions        = new COption*[iNumOfOptions];
    m_ppCClArgs         = new CClArg*[iNumOfOptions];
    m_iNumOfOptions     = iNumOfOptions;
    
    for (int i = 0; i < m_iNumOfOptions; i++)
    {
        m_ppCOptions[i] = new COption(myOptions[i]);
        m_ppCClArgs[i]  = new CClArg();
    }
 
    m_bIsInitialized    = true;

    return rErr;
}

Error_t CCommandLineOptions::resetInstance()
{
    
    m_bIsInitialized    = false;

    if (m_ppCOptions)
    {
        for (int i = 0; i < m_iNumOfOptions; i++)
            delete m_ppCOptions[i];
    }
    delete [] m_ppCOptions;
    m_ppCOptions        = 0;
    if (m_ppCClArgs)
    {
        for (int i = 0; i < m_iNumOfOptions; i++)
            delete m_ppCClArgs[i];
    }
    delete [] m_ppCClArgs;
    m_ppCClArgs         = 0;
    m_iNumOfOptions     = 0;

    return kNoError;

}

Error_t CCommandLineOptions::process( int argc, char* argv[] )
{
    if (!m_bIsInitialized)
        return kNotInitializedError;

    if (!argv || argc < 2)
        return kFunctionInvalidArgsError;

    for (int i = 1; i < argc; i++)
    {
        int iIdx = getBufferIdx(argv[i]);

        if (iIdx < 0)
            continue;

        assert(iIdx <= m_iNumOfOptions);

        if (m_ppCOptions[iIdx]->getType() == COption::OptionType_t::kBool)
            m_ppCClArgs[iIdx]->setCurrOption();
        else
        {
            if (i < argc-1)
            {
                m_ppCClArgs[iIdx]->setCurrOption(argv[i+1]);
                if (getBufferIdx(argv[i+1]) >= 0)
                {
                    return kInvalidString;
                }
                i++;
            }
        }
    }

    return kNoError;
}

int CCommandLineOptions::getBufferIdx( std::string argv ) const
{
    for (int i= 0; i < m_iNumOfOptions; i++)
    {
        if ((m_ppCOptions[i]->getClArg()).compare(argv) == 0)
            return i;
    }
    return -1;
}

int CCommandLineOptions::getBufferIdx( int iIdentifier ) const
{
    for (int i= 0; i < m_iNumOfOptions; i++)
    {
        if (m_ppCOptions[i]->getIdentifier() == iIdentifier)
            return i;
    }
    return -1;
}

bool CCommandLineOptions::isOptionSet( int iIdentifier ) const
{
    int iIdx = getBufferIdx(iIdentifier);

    if (iIdx >+ m_iNumOfOptions || iIdx < 0)
        return false;
    return m_ppCClArgs[iIdx]->isSet();
}

bool CCommandLineOptions::isOptionSet( std::string argv ) const
{
    int iIdx = getBufferIdx(argv);

    if (iIdx >= m_iNumOfOptions || iIdx < 0)
        return false;
    return m_ppCClArgs[iIdx]->isSet();
}

Error_t CCommandLineOptions::getOption( int iIdentifier, std::string& sResult ) const
{
    if (!isOptionSet(iIdentifier))
        return kNotInitializedError;

    int iIdx = getBufferIdx(iIdentifier);

    sResult = m_ppCClArgs[iIdx]->getArg();

    return kNoError;
}

Error_t CCommandLineOptions::getOption( std::string argv, std::string& sResult ) const
{
    if (!isOptionSet(argv))
        return kNotInitializedError;

    int iIdx = getBufferIdx(argv);

    sResult = m_ppCClArgs[iIdx]->getArg();

    return kNoError;
}
