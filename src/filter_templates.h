#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>

/*!
Templated variants of selected filters.
*/

/*!
Filter abstract base class.
*/
template <typename T>
class FilterBaseT {
public:
    virtual ~FilterBaseT() = default;
    virtual T filter_virtual(const T& input) = 0;
};


/*!
Null filter.
*/
template <typename T>
class FilterNullT : public FilterBaseT<T> {
public:
    FilterNullT() = default;
public:
    void init(float k) { (void)k; }
    void reset() {}
    void set_to_passthrough() {}
    void set_cutoff_frequency(float cutoff_frequency_hz, float dt) { (void)cutoff_frequency_hz; (void)dt; }
    void set_cutoff_frequency_and_reset(float cutoff_frequency_hz, float dt) { (void)cutoff_frequency_hz; (void)dt; }

    T filter(const T& input) { return input; }
    T filter(const T& input, float dt) { (void)dt; return input; }
    virtual T filter_virtual(const T& input) override { return filter(input); }
};


/*!
First order power transfer filter
*/
template <typename T>
class PowerTransferFilter1T : public FilterBaseT<T> {
public:
    explicit PowerTransferFilter1T(float k) : _k(k) {}
    PowerTransferFilter1T() : PowerTransferFilter1T(1.0F) {}
    PowerTransferFilter1T(float cutoff_frequency_hz, float dt) : PowerTransferFilter1T(gain_from_frequency(cutoff_frequency_hz, dt)) {}
public:
    void init(float k) { _k = k; reset(); }
    void reset() { _state = {}; }
    void set_to_passthrough() { _k = 1.0F; reset(); }

    T filter(const T& input) {
        _state = _state + (input - _state)*_k; // equivalent to _state = _k*input + (1.0F - _k)*_state;
        return _state;
    }
    virtual T filter_virtual(const T& input) override { return filter(input); }

    void set_cutoff_frequency(float cutoff_frequency_hz, float dt) { _k = gain_from_frequency(cutoff_frequency_hz, dt); }
    void set_cutoff_frequency_and_reset(float cutoff_frequency_hz, float dt) { _k = gain_from_frequency(cutoff_frequency_hz, dt); reset(); }
    // Calculates filter gain based on delay (time constant of filter) - time it takes for filter response to reach 63.2% of a step input.
    static float gain_from_delay(float delay, float dt) {
        if (delay <= 0) { return 1.0F; } // gain of 1.0F means no filtering
        const float omega = dt/delay;
        return omega/(omega + 1.0F);
    }
    static float gain_from_frequency(float cutoff_frequency_hz, float dt) {
        const float omega = 2.0F*PI_F*cutoff_frequency_hz*dt;
        return omega/(omega + 1.0F);
    }
// for testing
    const T& get_state() const { return _state; }
protected:
    float _k;
    T _state {};
protected:
    static constexpr float PI_F = 3.14159265358979323846F;
};


/*!
Second order power transfer filter
*/
template <typename T>
class PowerTransferFilter2T : public FilterBaseT<T> {
public:
    explicit PowerTransferFilter2T(float k) : _k(k) {}
    PowerTransferFilter2T() : PowerTransferFilter2T(1.0F) {}
    PowerTransferFilter2T(float cutoff_frequency_hz, float dt) : PowerTransferFilter2T(gain_from_frequency(cutoff_frequency_hz, dt)) {}
public:
    void init(float k) { _k = k; reset(); }
    void reset() { _state[0] = {}; _state[1] = {}; }
    void set_to_passthrough() { _k = 1.0F; }

    T filter(const T& input) {
        _state[1] = _state[1] + (input - _state[1])*_k;
        _state[0] = _state[0] + (_state[1] - _state[0])*_k;
        return _state[0];
    }
    virtual T filter_virtual(const T& input) override { return filter(input); }

    void set_cutoff_frequency(float cutoff_frequency_hz, float dt) { _k = gain_from_frequency(cutoff_frequency_hz, dt); }
    void set_cutoff_frequency_and_reset(float cutoff_frequency_hz, float dt) { _k = gain_from_frequency(cutoff_frequency_hz, dt); reset(); }
    static float gain_from_delay(float delay, float dt) {
        return PowerTransferFilter1T<T>::gain_from_delay(delay*CUTOFF_CORRECTION, dt);
    }
    static float gain_from_frequency(float cutoff_frequency_hz, float dt) {
        // shift cutoffFrequency to satisfy -3dB cutoff condition
        return PowerTransferFilter1T<T>::gain_from_frequency(cutoff_frequency_hz*CUTOFF_CORRECTION, dt);
    }
// for testing
    const std::array<T, 2>& get_state() const { return _state; }
protected:
    // PowerTransferFilter<n> cutoff correction = 1/sqrt(2^(1/n) - 1)
    static constexpr float CUTOFF_CORRECTION = 1.553773974F;
    float _k;
    std::array<T, 2> _state {};
};


/*!
Third order power transfer filter
*/
template <typename T>
class PowerTransferFilter3T : public FilterBaseT<T> {
public:
    explicit PowerTransferFilter3T(float k) : _k(k) {}
    PowerTransferFilter3T() : PowerTransferFilter3T(1.0F) {}
    PowerTransferFilter3T(float cutoff_frequency_hz, float dt) : PowerTransferFilter3T(gain_from_frequency(cutoff_frequency_hz, dt)) {}
public:
    void init(float k) { _k = k; reset(); }
    void reset() { _state[0] = {}; _state[1] = {}; _state[2] = {}; }
    void set_to_passthrough() { _k = 1.0F; reset(); }

    T filter(const T& input) {
        _state[2] = _state[2] + (input - _state[2])*_k;
        _state[1] = _state[1] + (_state[2] - _state[1])*_k;
        _state[0] = _state[0] + (_state[1] - _state[0])*_k;
        return _state[0];
    }
    virtual T filter_virtual(const T& input) override { return filter(input); }

    void set_cutoff_frequency(float cutoff_frequency_hz, float dt) { _k = gain_from_frequency(cutoff_frequency_hz, dt); }
    void set_cutoff_frequency_and_reset(float cutoff_frequency_hz, float dt) { _k = gain_from_frequency(cutoff_frequency_hz, dt); reset(); }
    static float gain_from_delay(float delay, float dt) {
        return PowerTransferFilter1T<T>::gain_from_delay(delay*CUTOFF_CORRECTION, dt);
    }
    static float gain_from_frequency(float cutoff_frequency_hz, float dt) {
        // shift cutoffFrequency to satisfy -3dB cutoff condition
        return PowerTransferFilter1T<T>::gain_from_frequency(cutoff_frequency_hz*CUTOFF_CORRECTION, dt);
    }
// for testing
    const std::array<T, 3>& get_state() const { return _state; }
protected:
    // PowerTransferFilter<n> cutoff correction = 1/sqrt(2^(1/n) - 1)
    static constexpr float CUTOFF_CORRECTION = 1.961459177F;
    float _k;
    std::array<T, 3> _state {};
};


/*!
Biquad filter, see https://en.wikipedia.org/wiki/Digital_biquad_filter

Has additional `_weight` member data, which allows the filter to combine input and output using `filter_weighted()` function.
`_weight` is ignored when using `filter()` function.
*/
template <typename T>
class BiquadFilterT : public FilterBaseT<T> {
public:
    BiquadFilterT(float a1, float a2, float b0, float b1, float b2) :
        _weight(1.0F),
        _a1(a1), _a2(a2),
        _b0(b0), _b1(b1), _b2(b2)
        {}
    BiquadFilterT() : BiquadFilterT(0.0F, 0.0F, 1.0F, 0.0F, 0.0F) {}
    struct state_t {
        T x1;
        T x2;
        T y1;
        T y2;
    };
public:
    void set_weight(float weight) { _weight = weight; }
    float get_weight() const { return _weight; }
    void set_parameters(float a1, float a2, float b0, float b1, float b2, float weight) {
        _weight = weight;
        _a1 = a1;
        _a2 = a2;
        _b0 = b0;
        _b1 = b1;
        _b2 = b2;
    }
    void set_parameters(float a1, float a2, float b0, float b1, float b2) {
        set_parameters(a1, a2, b0, b1, b2, 1.0F);
    }
    //! Copy parameters from another Biquad filter
    void set_parameters(const BiquadFilterT& other) {
        _weight = other._weight;
        _a1 = other._a1;
        _a2 = other._a2;
        _b0 = other._b0;
        _b1 = other._b1;
        _b2 = other._b2;
    }

    void reset() { _state.x1 = {}; _state.x2 = {}; _state.y1 = {}; _state.y2 = {}; }
    void set_to_passthrough() { _b0 = 1.0F; _b1 = 0.0F; _b2 = 0.0F; _a1 = 0.0F; _a2 = 0.0F;  _weight = 1.0F; reset(); }

    T filter(const T& input) {
        const T output = input*_b0 + _state.x1*_b1 + _state.x2*_b2 - _state.y1*_a1 - _state.y2*_a2;
        _state.x2 = _state.x1;
        _state.x1 = input;
        _state.y2 = _state.y1;
        _state.y1 = output;
        return output;
    }
    virtual T filter_virtual(const T& input) override { return filter(input); }

    T filter_weighted(const T& input) {
        const T output = filter(input);
        // weight of 1.0 gives just output, weight of 0.0 gives just input
        return (output - input)*_weight + input;
    }

    void init_lowpass(float frequency_hz, float looptime_seconds, float q) {
        assert(q != 0.0F && "q cannot be zero");
        set_looptime(looptime_seconds);
        set_q(q);
        set_low_pass_frequency(frequency_hz);
        reset();
    }
    void init_notch(float frequency_hz, float looptime_seconds, float q) {
        assert(q != 0.0F && "q cannot be zero");
        set_looptime(looptime_seconds);
        set_q(q);
        set_notch_frequency(frequency_hz);
        reset();
    }

    float calculate_omega(float frequency) const { return frequency*_2_pi_looptime_seconds; }

    void set_low_pass_frequencyWeighted(float frequency_hz, float weight);
    void set_low_pass_frequency(float frequency_hz) { set_low_pass_frequencyWeighted(frequency_hz, 1.0F); }

    void set_notch_frequency_weighted(float frequency_hz, float weight); // assumes q already set
    void set_notch_frequency(float frequency_hz) {set_notch_frequency_weighted(frequency_hz, 1.0F); } // assumes q already set
    void set_notch_frequency_weighted(float sin_omega, float two_cos_omega, float weight);
    void set_notch_frequency(float center_frequency_hz, float lower_cutoff_frequency_hz) {
        set_q(calculate_q(center_frequency_hz, lower_cutoff_frequency_hz));
        set_notch_frequency(center_frequency_hz);
    }
    void set_notch_frequency(uint16_t center_frequency_hz, uint16_t lower_cutoff_frequency_hz) {
        set_notch_frequency(static_cast<float>(center_frequency_hz), static_cast<float>(lower_cutoff_frequency_hz));
    }

    static float calculate_q(float center_frequency_hz, float lower_cutoff_frequency_hz) {
        return center_frequency_hz*lower_cutoff_frequency_hz/(center_frequency_hz*center_frequency_hz - lower_cutoff_frequency_hz*lower_cutoff_frequency_hz);
    }
    void set_q(float centerFrequency, float lower_cutoff_frequency) { _2q_reciprocal = 1.0F/(2.0F*calculate_q(centerFrequency, lower_cutoff_frequency)); }
    void set_q(float q) { _2q_reciprocal = 1.0F /(2.0F*q); }
    float get_q() const { return (1.0F/_2q_reciprocal)/2.0F; }

    void set_looptime(float looptime_seconds) { _2_pi_looptime_seconds = 2.0F*PI_F*looptime_seconds; }
// for testing
    const state_t& get_state() const { return _state; }
protected:
    float _weight {1.0F}; //<! weight of 1.0 gives just output, weight of 0.0 gives just input
    float _a1;
    float _a2;
    float _b0;
    float _b1;
    float _b2;

    state_t _state {};

    float _2q_reciprocal {1.0F}; // store 1/(2*q), since that is what is used in set_notch_frequency calculations
    float _2_pi_looptime_seconds {0.0F}; // store 2*PI*looptime_seconds, since that is what is used in calculations
protected:
    static constexpr float PI_F = 3.14159265358979323846F;
};

/*!
Note: weight must be in range [0, 1].
*/
template <typename T>
inline void BiquadFilterT<T>::set_low_pass_frequencyWeighted(float frequency_hz, float weight)
{
    _weight = weight;

    const float omega = frequency_hz*_2_pi_looptime_seconds;
#if defined(LIBRARY_FILTER_USE_SINCOS)
    float sin_omega {};
    float cos_omega {};
    sincosf(omega, &sin_omega, &cos_omega);
    const float alpha = sin_omega*_2q_reciprocal;
#else
    const float cos_omega = cosf(omega);
    const float alpha = sinf(omega)*_2q_reciprocal;
#endif
    const float a0_reciprocal = 1.0F/(1.0F + alpha);

    _b1 = (1.0F - cos_omega)*a0_reciprocal;
    _b0 = _b1*0.5F;
    _b2 = _b0;
    _a1 = -2.0F*cos_omega*a0_reciprocal;
    _a2 = (1.0F - alpha)*a0_reciprocal;
}

/*!
Note: weight must be in range [0, 1].
*/
template <typename T>
inline void BiquadFilterT<T>::set_notch_frequency_weighted(float frequency_hz, float weight)
{
    _weight = weight;

    const float omega = frequency_hz*_2_pi_looptime_seconds;
#if defined(LIBRARY_FILTER_USE_SINCOS)
    float sin_omega {};
    float cos_omega {};
    sincosf(omega, &sin_omega, &cos_omega);
    const float alpha = sin_omega*_2q_reciprocal;
#else
    const float cos_omega = cosf(omega);
    const float alpha = sinf(omega)*_2q_reciprocal;
#endif
    const float a0_reciprocal = 1.0F/(1.0F + alpha);

    _b0 = a0_reciprocal;
    _b2 = a0_reciprocal;
    _b1 = -2.0F*cos_omega*a0_reciprocal;
    _a1 = _b1;
    _a2 = (1.0F - alpha)*a0_reciprocal;
}

/*!
Note: weight must be in range [0, 1].
*/
template <typename T>
inline void BiquadFilterT<T>::set_notch_frequency_weighted(float sin_omega, float two_cos_omega, float weight)
{
    _weight = weight;

    const float alpha = sin_omega*_2q_reciprocal;
    const float a0_reciprocal = 1.0F/(1.0F + alpha);

    _b0 = a0_reciprocal;
    _b2 = a0_reciprocal;
    _b1 = -two_cos_omega*a0_reciprocal;
    _a1 = _b1;
    _a2 = (1.0F - alpha)*a0_reciprocal;
}


/*!
Simple moving average filter.
*/
template <typename T, size_t N>
class FilterMovingAverageT : public FilterBaseT<T> {
public:
    FilterMovingAverageT() {} // cppcheck-suppress uninitMemberVar
public:
    void reset() { _sum = {}; _count = 0; _index = 0;}

    T filter(const T& input);
    T filter(const T& input, float dt) { (void)dt; return filter(input); }
    virtual T filter_virtual(const T& input) override { return filter(input); }
protected:
    size_t _count {0};
    size_t _index {0};
    T _sum {};
    T _samples[N];
};

template <typename T, size_t N>
inline T FilterMovingAverageT<T, N>::filter(const T& input)
{
    _sum = _sum + input;
    if (_count < N) {
        _samples[_index] = input;
        ++_index;
        ++_count;
        return _sum*(1.0F/static_cast<float>(_count));
    }
    if (_index == N) {
        _index = 0;
    }
    _sum = _sum - _samples[_index];
    _samples[_index] = input;
    ++_index;
    static constexpr float N_RECIPROCAL = 1.0F/static_cast<float>(N);
    return _sum*N_RECIPROCAL;
}
