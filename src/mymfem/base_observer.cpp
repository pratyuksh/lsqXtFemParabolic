#include "base_observer.hpp"

#include <iostream>

using namespace mfem;


mymfem::BaseObserver
:: BaseObserver (const nlohmann::json& config)
{
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "visualization",
                                        m_boolVisualize, false);

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "dump_output",
                                        m_boolDumpOut, false);
}

mymfem::BaseObserver
:: BaseObserver (const nlohmann::json& config, int spatialLevel)
    : BaseObserver(config)
{
    m_meshNameSuffix = "_lx"+std::to_string(spatialLevel);
}

//! Uses GLVis, refer to examples in MFEM documentation.
void mymfem::BaseObserver
:: visualize (GridFunction& u) const
{
    socketstream sout;
    if (m_boolVisualize)
    {
        char vishost[] = "localhost";
        int  visport   = 19916;
        sout.open(vishost, visport);
        if (!sout.is_open())
        {
            std::cout << "Unable to connect to GLVis server at "
                      << vishost << ':' << visport << std::endl;
            m_boolVisualize = false;
            std::cout << "GLVis visualization disabled.\n";
        }
        else
        {
            sout.precision(m_precision);
            sout << "solution\n"
                 << *(u.FESpace()->GetMesh()) << u;
            sout << "pause\n";
            sout << std::flush;
        }
        sout.close();
    }
}

void mymfem::BaseObserver
:: dumpMesh (std::shared_ptr<Mesh>& mesh) const
{
    if (m_boolDumpOut)
    {
        std::string meshName
                = m_outputDir+"mesh"+m_meshNameSuffix;
        std::cout << meshName << std::endl;

        std::ofstream meshOfs(meshName.c_str());
        meshOfs.precision(m_precision);
        mesh->Print(meshOfs);
        meshOfs.close();
    }
}


// End of file
