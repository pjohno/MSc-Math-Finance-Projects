#pragma once
#include "xeus/xjson.hpp"

namespace im
{
    xeus::xjson mime_bundle_repr(const MSC_PROJECTS::gnuplotImage& i)
    {
        auto bundle = xeus::xjson::object();
        bundle["image/png"] = xtl::base64encode(i.imageText);
        return bundle;
    }
}