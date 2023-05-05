/*
 * MOD Convolution Loader
 * Copyright (C) 2022-2023 Filipe Coelho <falktx@falktx.com>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

// IDE helper (not needed for building)
#include "DistrhoPluginInfo.h"

#include "DistrhoPlugin.hpp"
#include "Korg35Filters.hpp"

#include "extra/ScopedPointer.hpp"
#include "extra/ValueSmoother.hpp"

#include "dr_flac.h"
#include "dr_wav.h"
// -Wunused-variable
#include "r8brain/CDSPResampler.h"

#include "TwoStageThreadedConvolver.hpp"

START_NAMESPACE_DISTRHO

#if defined(_MOD_DEVICE_DUO)
static constexpr const size_t headBlockSize = 256;
static constexpr const size_t tailBlockSize = 4096;
#elif defined(_MOD_DEVICE_DWARF)
static constexpr const size_t headBlockSize = 128;
static constexpr const size_t tailBlockSize = 2048;
#else
static constexpr const size_t headBlockSize = 128;
static constexpr const size_t tailBlockSize = 1024;
#endif

// -----------------------------------------------------------------------

class OneKnobConvolutionReverbPlugin : public Plugin
{
public:
    OneKnobConvolutionReverbPlugin()
        : Plugin(kParameterCount, kProgramCount, kStateCount)
    {
        for (uint i=0; i<kParameterCount; ++i)
            parameters[i] = kParameterRanges[i].def;

        const float sampleRate = static_cast<float>(getSampleRate());

        korgFilterL.setSampleRate(sampleRate);
        korgFilterR.setSampleRate(sampleRate);

        korgFilterL.setFrequency(kParameterRanges[kParameterHighPassFilter].def);
        korgFilterR.setFrequency(kParameterRanges[kParameterHighPassFilter].def);

        smoothDryLevel.setSampleRate(sampleRate);
        smoothWetLevel.setSampleRate(sampleRate);

        smoothDryLevel.setTimeConstant(0.1f);
        smoothWetLevel.setTimeConstant(0.1f);

        smoothDryLevel.setTargetValue(std::pow(10.f, 0.05f * kParameterRanges[kParameterDryLevel].def));
        smoothWetLevel.setTargetValue(std::pow(10.f, 0.05f * kParameterRanges[kParameterWetLevel].def));
    }

    ~OneKnobConvolutionReverbPlugin() override
    {
    }

protected:
    // -------------------------------------------------------------------
    // Information

    const char* getLabel() const noexcept override
    {
        return DISTRHO_PLUGIN_NAME;
    }

    const char* getMaker() const noexcept override
    {
        return DISTRHO_PLUGIN_BRAND;
    }

    const char* getHomePage() const override
    {
        return DISTRHO_PLUGIN_URI;
    }

    uint32_t getVersion() const noexcept override
    {
        return d_version(1, 0, 0);
    }

    const char* getDescription() const override
    {
        return "The MOD Convolution Loader enables you to easily create custom reverb and other effects using impulse responses (IRs)";
    }

    const char* getLicense() const noexcept override
    {
        return "ISC";
    }

    int64_t getUniqueId() const noexcept override
    {
        return d_cconst('O', 'K', 'c', 'r');
    }

    // -------------------------------------------------------------------
    // Init

    void initParameter(uint32_t index, Parameter& parameter) override
    {
        switch (index)
        {
        case kParameterDryLevel:
            parameter.hints = kParameterIsAutomatable;
            parameter.name = "Dry Level";
            parameter.symbol = "drylevel";
            parameter.unit = "dB";
            parameter.ranges.def = kParameterRanges[kParameterDryLevel].def;
            parameter.ranges.min = kParameterRanges[kParameterDryLevel].min;
            parameter.ranges.max = kParameterRanges[kParameterDryLevel].max;
            {
                ParameterEnumerationValue* const enumValues =  new ParameterEnumerationValue[1];
                enumValues[0].value = kParameterRanges[kParameterDryLevel].min;
                enumValues[0].label = "Off";
                parameter.enumValues.count = 1;
                parameter.enumValues.values = enumValues;
            }
            break;
        case kParameterWetLevel:
            parameter.hints = kParameterIsAutomatable;
            parameter.name = "Wet Level";
            parameter.symbol = "wetlevel";
            parameter.unit = "dB";
            parameter.ranges.def = kParameterRanges[kParameterWetLevel].def;
            parameter.ranges.min = kParameterRanges[kParameterWetLevel].min;
            parameter.ranges.max = kParameterRanges[kParameterWetLevel].max;
            {
                ParameterEnumerationValue* const enumValues = new ParameterEnumerationValue[1];
                enumValues[0].value = kParameterRanges[kParameterWetLevel].min;
                enumValues[0].label = "Off";
                parameter.enumValues.count = 1;
                parameter.enumValues.values = enumValues;
            }
            break;
        case kParameterHighPassFilter:
            parameter.hints = kParameterIsAutomatable;
            parameter.name = "High Pass Filter";
            parameter.symbol = "hpf";
            parameter.unit = "Hz";
            parameter.ranges.def = kParameterRanges[kParameterHighPassFilter].def;
            parameter.ranges.min = kParameterRanges[kParameterHighPassFilter].min;
            parameter.ranges.max = kParameterRanges[kParameterHighPassFilter].max;
            {
                ParameterEnumerationValue* const enumValues = new ParameterEnumerationValue[1];
                enumValues[0].value = 0.f;
                enumValues[0].label = "Off";
                parameter.enumValues.count = 1;
                parameter.enumValues.values = enumValues;
            }
            break;
        case kParameterTrails:
            parameter.hints = kParameterIsAutomatable | kParameterIsInteger | kParameterIsBoolean;
            parameter.name = "Trails";
            parameter.symbol = "trails";
            parameter.ranges.def = kParameterRanges[kParameterTrails].def;
            parameter.ranges.min = kParameterRanges[kParameterTrails].min;
            parameter.ranges.max = kParameterRanges[kParameterTrails].max;
            break;
        case kParameterBypass:
            parameter.initDesignation(kParameterDesignationBypass);
            break;
        case kParameterBuffered:
            parameter.hints = kParameterIsAutomatable | kParameterIsInteger | kParameterIsBoolean;
            parameter.name = "Buffered";
            parameter.symbol = "buffered";
            parameter.ranges.def = kParameterRanges[kParameterTrails].def;
            parameter.ranges.min = kParameterRanges[kParameterTrails].min;
            parameter.ranges.max = kParameterRanges[kParameterTrails].max;
            break;
        }
    }

    void initState(uint32_t index, State &state) override
    {
        switch (index)
        {
        case kStateFile:
            state.hints = kStateIsFilenamePath;
            state.key = "irfile";
            state.label = "IR File";
           #ifdef __MOD_DEVICES__
            state.fileTypes = "ir";
           #endif
            break;
        }
    }

    // -------------------------------------------------------------------
    // Parameters

    float getParameterValue(const uint32_t index) const override
    {
        return parameters[index];
    }

    void setParameterValue(const uint32_t index, const float value) override
    {
        parameters[index] = value;

        switch (index)
        {
        case kParameterDryLevel:
            if (!bypassed)
                smoothDryLevel.setTargetValue(std::pow(10.f, 0.05f * value));
            break;
        case kParameterWetLevel:
            if (!bypassed)
                smoothWetLevel.setTargetValue(std::pow(10.f, 0.05f * value));
            break;
        case kParameterHighPassFilter:
            korgFilterL.setFrequency(value);
            korgFilterR.setFrequency(value);
            break;
        case kParameterTrails:
            trails = value > 0.5f;
            if (bypassed)
                smoothWetLevel.setTargetValue(trails ? std::pow(10.f, 0.05f * parameters[kParameterWetLevel]) : 0.f);
            break;
        case kParameterBypass:
            bypassed = value > 0.5f;
            if (bypassed)
            {
                smoothDryLevel.setTargetValue(1.f);
                smoothWetLevel.setTargetValue(trails ? std::pow(10.f, 0.05f * parameters[kParameterWetLevel]) : 0.f);
            }
            else
            {
                korgFilterL.reset();
                korgFilterR.reset();
                smoothDryLevel.setTargetValue(std::pow(10.f, 0.05f * parameters[kParameterDryLevel]));
                smoothWetLevel.setTargetValue(std::pow(10.f, 0.05f * parameters[kParameterWetLevel]));
            }
            break;
        case kParameterBuffered:
            break;
        }
    }

    void setState(const char* const key, const char* const value) override
    {
        if (std::strcmp(key, "irfile") == 0)
        {
            unsigned int channels;
            unsigned int sampleRate;
            drwav_uint64 numFrames;
            const size_t valuelen = std::strlen(value);

            ScopedPointer<TwoStageThreadedConvolver> newConvolverL, newConvolverR;

            if (valuelen <= 5)
            {
                bufferedConvolver.stop();
                const MutexLocker cml(mutex);
                convolverL.swapWith(newConvolverL);
                convolverR.swapWith(newConvolverR);
                return;
            }

            float* ir;
            if (::strncasecmp(value + (std::max(size_t(0), valuelen - 5u)), ".flac", 5) == 0)
                ir = drflac_open_file_and_read_pcm_frames_f32(value, &channels, &sampleRate, &numFrames, nullptr);
            else
                ir = drwav_open_file_and_read_pcm_frames_f32(value, &channels, &sampleRate, &numFrames, nullptr);
            DISTRHO_SAFE_ASSERT_RETURN(ir != nullptr,);

            loadedFilename = value;

            float* irBufL;
            float* irBufR;
            switch (channels)
            {
            case 1:
                irBufL = irBufR = ir;
                break;
            case 2:
                irBufL = new float[numFrames];
                irBufR = new float[numFrames];
                for (drwav_uint64 i = 0, j = 0; i < numFrames; ++i)
                {
                    irBufL[i] = ir[j++];
                    irBufR[i] = ir[j++];
                }
                break;
            case 4:
                irBufL = new float[numFrames];
                irBufR = new float[numFrames];
                for (drwav_uint64 i = 0, j = 0; i < numFrames; ++i, j += 4)
                {
                    irBufL[i] = ir[j + 0] + ir[j + 2];
                    irBufR[i] = ir[j + 1] + ir[j + 3];
                }
                break;
            default:
                irBufL = new float[numFrames];
                irBufR = new float[numFrames];
                for (drwav_uint64 i = 0, j = 0; i < numFrames; ++i)
                {
                    irBufL[i] = irBufR[i] = ir[j];
                    j += channels;
                }
                break;
            }

            if (sampleRate != getSampleRate())
            {
                r8b::CDSPResampler16IR resampler(sampleRate, getSampleRate(), numFrames);
                const int numResampledFrames = resampler.getMaxOutLen(0);
                DISTRHO_SAFE_ASSERT_RETURN(numResampledFrames > 0,);

                // left channel, always present
                float* const irBufResampledL = new float[numResampledFrames];
                resampler.oneshot(irBufL, numFrames, irBufResampledL, numResampledFrames);
                delete[] irBufL;
                irBufL = irBufResampledL;

                // right channel, optional
                if (irBufL != irBufR)
                {
                    float* const irBufResampledR = new float[numResampledFrames];
                    resampler.oneshot(irBufR, numFrames, irBufResampledR, numResampledFrames);
                    delete[] irBufR;
                    irBufR = irBufResampledR;
                }

                numFrames = numResampledFrames;
            }

            newConvolverL = new TwoStageThreadedConvolver();
            newConvolverL->init(headBlockSize, tailBlockSize, irBufL, numFrames);

            newConvolverR = new TwoStageThreadedConvolver();
            newConvolverR->init(headBlockSize, tailBlockSize, irBufR, numFrames);

            bufferedConvolver.stop();
            bufferedConvolver.start(newConvolverL, newConvolverR);

            {
                const MutexLocker cml(mutex);
                convolverL.swapWith(newConvolverL);
                convolverR.swapWith(newConvolverR);
            }

            if (irBufL != ir)
                delete[] irBufL;
            if (irBufR != irBufL)
                delete[] irBufR;

            drwav_free(ir, nullptr);
            return;
        }
    }

    // -------------------------------------------------------------------
    // Process

    void activate() override
    {
        const uint32_t bufSize = bufferSize = getBufferSize();

        highpassBufL = new float[bufSize];
        highpassBufR = new float[bufSize];
        inplaceProcBufL = new float[bufSize];
        inplaceProcBufR = new float[bufSize];

        korgFilterL.reset();
        korgFilterR.reset();

        smoothDryLevel.clearToTargetValue();
        smoothWetLevel.clearToTargetValue();
    }

    void deactivate() override
    {
        delete[] highpassBufL;
        delete[] highpassBufR;
        delete[] inplaceProcBufL;
        delete[] inplaceProcBufR;
        bufferSize = 0;
        highpassBufL = highpassBufR = nullptr;
        inplaceProcBufL = inplaceProcBufR = nullptr;
    }

    void run(const float** const inputs, float** const outputs, const uint32_t frames) override
    {
        if (frames == 0)
            return;

        const float* const inL = inputs[0];
        const float* const inR = inputs[1];
        /* */ float* const outL = outputs[0];
        /* */ float* const outR = outputs[1];

        // optimize for non-denormal usage
        for (uint32_t i = 0; i < frames; ++i)
        {
            if (!std::isfinite(inL[i]))
                __builtin_unreachable();
            if (!std::isfinite(inR[i]))
                __builtin_unreachable();
            if (!std::isfinite(outL[i]))
                __builtin_unreachable();
            if (!std::isfinite(outR[i]))
                __builtin_unreachable();
        }

        const float* dryBufL = inL;
        const float* dryBufR = inR;

        const int hpf = static_cast<int>(parameters[kParameterHighPassFilter] + 0.5f);

        if (bypassed)
        {
            std::memset(highpassBufL, 0, sizeof(float) * frames);
            std::memset(highpassBufR, 0, sizeof(float) * frames);
        }
        else if (hpf == 0)
        {
            std::memcpy(highpassBufL, inL, sizeof(float) * frames);
            std::memcpy(highpassBufR, inR, sizeof(float) * frames);
        }
        else
        {
            korgFilterL.processHighPass(inL, highpassBufL, frames);
            korgFilterR.processHighPass(inR, highpassBufR, frames);
        }

        if (outL == inL)
        {
            dryBufL = inplaceProcBufL;
            std::memcpy(inplaceProcBufL, inL, sizeof(float) * frames);
        }

        if (outR == inR)
        {
            dryBufR = inplaceProcBufR;
            std::memcpy(inplaceProcBufR, inR, sizeof(float) * frames);
        }

        float wetLevel, dryLevel;

        const MutexTryLocker cmtl(mutex);

        if (cmtl.wasLocked())
        {
            TwoStageThreadedConvolver* const convL = convolverL.get();
            TwoStageThreadedConvolver* const convR = convolverR.get();

            if (convL != nullptr && convR != nullptr)
            {
                if (1)
                {
                    const float* const ins[2] = { highpassBufL, highpassBufR };
                    float* const outs[2] = { outL, outR };
                    bufferedConvolver.process(ins, outs, frames);
                }
                else
                {
                    convL->process(highpassBufL, outL, frames);
                    convR->process(highpassBufR, outR, frames);
                }

                for (uint32_t i = 0; i < frames; ++i)
                {
                    dryLevel = smoothDryLevel.next();
                    wetLevel = smoothWetLevel.next();

                    if (wetLevel <= 0.001f)
                    {
                        outL[i] = outR[i] = 0.f;
                    }
                    else
                    {
                        outL[i] *= wetLevel;
                        outR[i] *= wetLevel;
                    }

                    if (dryLevel > 0.001f)
                    {
                        outL[i] += dryBufL[i] * dryLevel;
                        outR[i] += dryBufR[i] * dryLevel;
                    }
                }

                return;
            }
        }

        for (uint32_t i = 0; i < frames; ++i)
        {
            smoothWetLevel.next();
            dryLevel = smoothDryLevel.next();

            outL[i] = dryBufL[i] * dryLevel;
            outR[i] = dryBufR[i] * dryLevel;
        }
    }

    void sampleRateChanged(const double newSampleRate) override
    {
        korgFilterL.setSampleRate(newSampleRate);
        korgFilterR.setSampleRate(newSampleRate);

        smoothDryLevel.setSampleRate(newSampleRate);
        smoothWetLevel.setSampleRate(newSampleRate);

        // reload file
        if (char* const filename = loadedFilename.getAndReleaseBuffer())
        {
            setState("irfile", filename);
            std::free(filename);
        }
    }

  // -------------------------------------------------------------------

private:
    ScopedPointer<TwoStageThreadedConvolver> convolverL, convolverR;
    StereoBufferedConvolver bufferedConvolver;
    Korg35Filter korgFilterL, korgFilterR;
    Mutex mutex;
    String loadedFilename;

    bool bypassed = false;
    bool trails = true;
    uint32_t bufferSize = 0;

    float parameters[kParameterCount];

    // smoothed parameters
    LinearValueSmoother smoothDryLevel;
    LinearValueSmoother smoothWetLevel;

    // buffers for placing highpass signal before convolution
    float* highpassBufL = nullptr;
    float* highpassBufR = nullptr;

    // if doing inline processing, copy buffers here before convolution
    float* inplaceProcBufL = nullptr;
    float* inplaceProcBufR = nullptr;

    DISTRHO_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(OneKnobConvolutionReverbPlugin)
};

// -----------------------------------------------------------------------

Plugin* createPlugin()
{
    return new OneKnobConvolutionReverbPlugin();
}

END_NAMESPACE_DISTRHO
