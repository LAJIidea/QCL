//
// Created by BY210033 on 2022/8/25.
//

#ifndef QCL_RANDOM_ENGINE_H
#define QCL_RANDOM_ENGINE_H

#include <random>

namespace qram_engine {

    struct random_engine
    {
        unsigned int seed = 10101;
        std::default_random_engine reng;
        random_engine() {}

        inline static random_engine& get_instance()
        {
            static random_engine inst;
            return inst;
        }

        inline std::default_random_engine& _get_engine()
        {
            return reng;
        }

        inline static std::default_random_engine& get_engine()
        {
            return get_instance()._get_engine();
        }

        inline double _rng()
        {
            static std::uniform_real_distribution<double> ud(0, 1);
            return ud(reng);
        }

        inline double uniform01()
        {
            return _rng();
        }

        inline void set_seed(unsigned int _seed)
        {
            seed = _seed;
            reng.seed(_seed);
        }

        inline unsigned int reseed()
        {
            set_seed((unsigned int)(_rng() * std::numeric_limits<unsigned int>::max()));
            return seed;
        }

        inline unsigned int get_seed() {
            return seed;
        }

        inline static unsigned int time_seed()
        {
            unsigned int seed = time(0);
            get_instance().set_seed(seed);
            return seed;
        }

    };

} // namespace qram_simulator

#endif //QCL_RANDOM_ENGINE_H
