#pragma once

#include <unordered_set>
#include <memory>

struct fluid_class
{
    std::unordered_set<std::shared_ptr<fluid_class>> soluble_classes;
};