/**
 * @author Tim Luchterhand
 * @date 06.05.2024
 * @brief CRTP helper class
 */

#ifndef TEMPO_CRTP_HPP
#define TEMPO_CRTP_HPP

namespace tempo {
    /**
     * @brief Helper class for CRTP base classes.
     * @details @copybrief
     * Provides accessors to the derived class.
     * @tparam Implementation Concrete implementation of the interface
     * @tparam CRTPImpl CRTP base class implementation
     * @tparam CRTPImplTArgs template arguments for the CRTP base class implementation
     */
    template<typename Implementation, template<typename...> typename CRTPImpl, typename ...CRTPImplTArgs>
    struct crtp {
        /**
         * Mutable access to concrete implementation
         * @return
         */
        Implementation &getImpl() { return static_cast<Implementation&>(*this); }

        /**
         * Const access to concrete implementation
         * @return
         */
        const Implementation &getImpl() const { return static_cast<const Implementation&>(*this); }

    private:
        crtp() = default;
        friend CRTPImpl<Implementation, CRTPImplTArgs...>;
    };
}

#endif //TEMPO_CRTP_HPP
