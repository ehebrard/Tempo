/**
 * @author Tim Luchterhand
 * @date 06.05.2024
 * @brief
 */

#ifndef TEMPO_CRTP_HPP
#define TEMPO_CRTP_HPP

namespace tempo {
    template<typename Implementation, template<typename...> typename CRTPImpl, typename ...ImplTArgs>
    struct crtp {
        Implementation &getImpl() { return static_cast<Implementation&>(*this); }

        const Implementation &getImpl() const { return static_cast<const Implementation&>(*this); }

    private:
        crtp() = default;
        friend CRTPImpl<Implementation, ImplTArgs...>;
    };
}

#endif //TEMPO_CRTP_HPP
