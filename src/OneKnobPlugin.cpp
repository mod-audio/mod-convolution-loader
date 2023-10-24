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

#if defined(_MOD_DEVICE_DUO) && defined(CONVOLUTION_REVERB)
static constexpr const size_t headBlockSize = 256;
static constexpr const size_t tailBlockSize = 4096;
#elif defined(_MOD_DEVICE_DWARF) && defined(CONVOLUTION_REVERB)
static constexpr const size_t headBlockSize = 128;
static constexpr const size_t tailBlockSize = 2048;
#else
static constexpr const size_t headBlockSize = 128;
static constexpr const size_t tailBlockSize = 1024;
#endif

// -----------------------------------------------------------------------

class ConvolutionLoaderPlugin : public Plugin
{
public:
    ConvolutionLoaderPlugin()
        : Plugin(kParameterCount, kProgramCount, kStateCount)
    {
        for (uint i=0; i<kParameterCount; ++i)
            parameters[i] = kParameterRanges[i].def;

       #ifdef CONVOLUTION_REVERB
        korgFilterL.setFrequency(kParameterRanges[kParameterHighPassFilter].def);
        korgFilterR.setFrequency(kParameterRanges[kParameterHighPassFilter].def);
       #endif

        smoothDryLevel.setTimeConstant(0.1f);
       #ifdef CONVOLUTION_REVERB
        smoothDryLevel.setTargetValue(std::pow(10.f, 0.05f * kParameterRanges[kParameterDryLevel].def));
       #else
        smoothDryLevel.setTargetValue(0.f);
       #endif

        smoothWetLevel.setTimeConstant(0.1f);
        smoothWetLevel.setTargetValue(std::pow(10.f, 0.05f * kParameterRanges[kParameterWetLevel].def));

        // init buffers
        bufferSizeChanged(getBufferSize());
        sampleRateChanged(getSampleRate());
    }

    ~ConvolutionLoaderPlugin() override
    {
        delete[] highpassBufL;
        delete[] inplaceProcBufL;
       #ifdef CONVOLUTION_REVERB
        delete[] highpassBufR;
        delete[] inplaceProcBufR;
       #endif
    }

protected:
    // -------------------------------------------------------------------
    // Information

    const char* getLabel() const noexcept override
    {
        return DISTRHO_PLUGIN_LABEL;
    }

    const char* getMaker() const noexcept override
    {
        return DISTRHO_PLUGIN_BRAND;
    }

    const char* getHomePage() const override
    {
        return "https://mod.audio/";
    }

    uint32_t getVersion() const noexcept override
    {
        return d_version(1, 0, 0);
    }

    const char* getDescription() const override
    {
       #ifdef CONVOLUTION_REVERB
        return R"(
The MOD Convolution Loader enables you to easily create custom reverb and other effects using impulse responses (IRs).

Using a very resource effective processing, it allows IRs with a length of up to 15s on the MOD Dwarf and MOD Duo X, less on the MOD Duo.

IRs can be uploaded using the file manager in either WAV and FLAC format.
Multi-channel files are supported but always processed in stereo mode.
For faster loading it is recommended to use WAV format files with 48kHz and 24 bit.

The plugin is very easy to use: Simply select the IR file, set the levels as desired, optionally filter out some low end with the built-in high pass filter, and off you go.

There is a Trails switch in the advanced settings;
when enabled the reverb tail peacefully decays when you bypass the plugin, disable it and the tail will immediately disappear when bypassed.

Features:
Plugin by MOD Audio & DISTRHO, based on HiFi-LoFi FFTConvolver engine.
)";
       #else
        return R"(
The MOD Cabinet Loader enables you to play through your favorite cabinet impulse responses (IRs).

The plugin is very easy to use: Simply select the IR file, set the level as desired and off you go.

IRs can be uploaded using the file manager in either WAV and FLAC format.
Multi-channel files are supported but always processed in mono mode using the first channel.
For faster loading it is recommended to use WAV format files with 48kHz and 24 bit.

Features:
Plugin by MOD Audio & DISTRHO, based on HiFi-LoFi FFTConvolver engine.
)";
       #endif
    }

    const char* getLicense() const noexcept override
    {
        return "ISC";
    }

    int64_t getUniqueId() const noexcept override
    {
       #ifdef CONVOLUTION_REVERB
        return d_cconst('M', 'A', 'c', 'r');
       #else
        return d_cconst('M', 'A', 'c', 'a');
       #endif
    }

    // -------------------------------------------------------------------
    // Init

    void initParameter(uint32_t index, Parameter& parameter) override
    {
        switch (index)
        {
       #ifdef CONVOLUTION_REVERB
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
       #endif
        case kParameterWetLevel:
            parameter.hints = kParameterIsAutomatable;
           #ifdef CONVOLUTION_REVERB
            parameter.name = "Wet Level";
            parameter.symbol = "wetlevel";
           #else
            parameter.name = "Output";
            parameter.symbol = "level";
           #endif
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
       #ifdef CONVOLUTION_REVERB
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
       #endif
        case kParameterBypass:
            parameter.initDesignation(kParameterDesignationBypass);
            break;
       #ifdef CONVOLUTION_REVERB
        case kParameterBuffered:
            parameter.hints = kParameterIsAutomatable | kParameterIsInteger | kParameterIsBoolean;
            parameter.name = "Buffered";
            parameter.symbol = "buffered";
            parameter.ranges.def = kParameterRanges[kParameterBuffered].def;
            parameter.ranges.min = kParameterRanges[kParameterBuffered].min;
            parameter.ranges.max = kParameterRanges[kParameterBuffered].max;
            break;
       #endif
        }
    }

    void initState(uint32_t index, State& state) override
    {
        switch (index)
        {
        case kStateFile:
            state.hints = kStateIsFilenamePath;
            state.key = "irfile";
            state.label = "IR File";
          #ifdef __MOD_DEVICES__
           #ifdef CONVOLUTION_REVERB
            state.fileTypes = "ir";
           #else
            state.fileTypes = "cabsim";
           #endif
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
       #ifdef CONVOLUTION_REVERB
        case kParameterDryLevel:
            if (!bypassed)
                smoothDryLevel.setTargetValue(std::pow(10.f, 0.05f * value));
            break;
       #endif
        case kParameterWetLevel:
            if (!bypassed)
                smoothWetLevel.setTargetValue(std::pow(10.f, 0.05f * value));
            break;
       #ifdef CONVOLUTION_REVERB
        case kParameterHighPassFilter:
            korgFilterL.setFrequency(value);
            korgFilterR.setFrequency(value);
            break;
        case kParameterTrails:
            trails = value > 0.5f;
            if (bypassed)
                smoothWetLevel.setTargetValue(trails ? std::pow(10.f, 0.05f * parameters[kParameterWetLevel]) : 0.f);
            break;
       #endif
        case kParameterBypass:
            bypassed = value > 0.5f;
            if (bypassed)
            {
                smoothDryLevel.setTargetValue(1.f);
                smoothWetLevel.setTargetValue(trails ? std::pow(10.f, 0.05f * parameters[kParameterWetLevel]) : 0.f);
            }
            else
            {
               #ifdef CONVOLUTION_REVERB
                korgFilterL.reset();
                korgFilterR.reset();
                smoothDryLevel.setTargetValue(std::pow(10.f, 0.05f * parameters[kParameterDryLevel]));
               #else
                smoothDryLevel.setTargetValue(0.f);
               #endif
                smoothWetLevel.setTargetValue(std::pow(10.f, 0.05f * parameters[kParameterWetLevel]));
            }
            break;
       #ifdef CONVOLUTION_REVERB
        case kParameterBuffered:
            buffered = value > 0.5f;
            break;
       #endif
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

           #ifdef CONVOLUTION_REVERB
            ScopedPointer<TwoStageThreadedConvolver> newConvolverL;
            ScopedPointer<TwoStageThreadedConvolver> newConvolverR;
           #else
            ScopedPointer<TwoStageThreadedConvolver> newConvolver;
           #endif

            if (valuelen <= 5)
            {
               #ifdef CONVOLUTION_REVERB
                bufferedConvolver.stop();
                const MutexLocker cml(mutex);
                convolverL.swapWith(newConvolverL);
                convolverR.swapWith(newConvolverR);
               #else
                const MutexLocker cml(mutex);
                convolver.swapWith(newConvolver);
               #endif
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
           #ifdef CONVOLUTION_REVERB
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
           #else
            if (channels == 1)
            {
                irBufL = ir;
            }
            else
            {
                irBufL = new float[numFrames];
                for (drwav_uint64 i = 0, j = 0; i < numFrames; ++i)
                {
                    irBufL[i] = ir[j];
                    j += channels;
                }
            }
           #endif

            if (sampleRate != getSampleRate())
            {
                r8b::CDSPResampler16IR resampler(sampleRate, getSampleRate(), numFrames);
                const int numResampledFrames = resampler.getMaxOutLen(0);
                DISTRHO_SAFE_ASSERT_RETURN(numResampledFrames > 0,);

                // left channel, always present
                float* const irBufResampledL = new float[numResampledFrames];
                resampler.oneshot(irBufL, numFrames, irBufResampledL, numResampledFrames);
                if (irBufL != ir)
                    delete[] irBufL;
                irBufL = irBufResampledL;

               #ifdef CONVOLUTION_REVERB
                // right channel, optional
                if (irBufL != irBufR)
                {
                    float* const irBufResampledR = new float[numResampledFrames];
                    resampler.oneshot(irBufR, numFrames, irBufResampledR, numResampledFrames);
                    if (irBufR != ir)
                        delete[] irBufR;
                    irBufR = irBufResampledR;
                }
               #endif

                numFrames = numResampledFrames;
            }

           #ifdef CONVOLUTION_REVERB
            newConvolverL = new TwoStageThreadedConvolver();
            newConvolverL->init(headBlockSize, tailBlockSize, irBufL, numFrames);

            newConvolverR = new TwoStageThreadedConvolver();
            newConvolverR->init(headBlockSize, tailBlockSize, irBufR, numFrames);

            bufferedConvolver.stop();
            bufferedConvolver.start(newConvolverL, newConvolverR);
           #else
            newConvolver = new TwoStageThreadedConvolver();
            newConvolver->init(headBlockSize, tailBlockSize, irBufL, numFrames);
           #endif

            {
                const MutexLocker cml(mutex);
               #ifdef CONVOLUTION_REVERB
                convolverL.swapWith(newConvolverL);
                convolverR.swapWith(newConvolverR);
               #else
                convolver.swapWith(newConvolver);
               #endif
            }

            if (irBufL != ir)
                delete[] irBufL;
           #ifdef CONVOLUTION_REVERB
            if (irBufR != irBufL)
                delete[] irBufR;
           #endif

            drwav_free(ir, nullptr);
            return;
        }
    }

    // -------------------------------------------------------------------
    // Process

    void activate() override
    {
       #ifdef CONVOLUTION_REVERB
        korgFilterL.reset();
        korgFilterR.reset();
       #endif

        smoothDryLevel.clearToTargetValue();
        smoothWetLevel.clearToTargetValue();
    }

    void run(const float** const inputs, float** const outputs, const uint32_t frames) override
    {
        if (frames == 0)
            return;

        const float* const inL  = inputs[0];
        /* */ float* const outL = outputs[0];
       #ifdef CONVOLUTION_REVERB
        const float* const inR  = inputs[1];
        /* */ float* const outR = outputs[1];
       #endif

        // optimize for non-denormal usage
        for (uint32_t i = 0; i < frames; ++i)
        {
            if (!std::isfinite(inL[i]))
                __builtin_unreachable();
            if (!std::isfinite(outL[i]))
                __builtin_unreachable();
           #ifdef CONVOLUTION_REVERB
            if (!std::isfinite(inR[i]))
                __builtin_unreachable();
            if (!std::isfinite(outR[i]))
                __builtin_unreachable();
           #endif
        }

        const float* dryBufL = inL;
       #ifdef CONVOLUTION_REVERB
        const float* dryBufR = inR;

        if (bypassed)
        {
            std::memset(highpassBufL, 0, sizeof(float) * frames);
            std::memset(highpassBufR, 0, sizeof(float) * frames);
        }
        else if (static_cast<int>(parameters[kParameterHighPassFilter] + 0.5f) == 0)
        {
            std::memcpy(highpassBufL, inL, sizeof(float) * frames);
            std::memcpy(highpassBufR, inR, sizeof(float) * frames);
        }
        else
        {
            korgFilterL.processHighPass(inL, highpassBufL, frames);
            korgFilterR.processHighPass(inR, highpassBufR, frames);
        }
       #else
        if (bypassed)
        {
            std::memset(highpassBufL, 0, sizeof(float) * frames);
        }
        else
        {
            std::memcpy(highpassBufL, inL, sizeof(float) * frames);
        }
       #endif

        if (outL == inL)
        {
            dryBufL = inplaceProcBufL;
            std::memcpy(inplaceProcBufL, inL, sizeof(float) * frames);
        }

       #ifdef CONVOLUTION_REVERB
        if (outR == inR)
        {
            dryBufR = inplaceProcBufR;
            std::memcpy(inplaceProcBufR, inR, sizeof(float) * frames);
        }
       #endif

        float wetLevel, dryLevel;

        const MutexTryLocker cmtl(mutex);

        if (cmtl.wasLocked())
        {
           #ifdef CONVOLUTION_REVERB
            TwoStageThreadedConvolver* const convL = convolverL.get();
            TwoStageThreadedConvolver* const convR = convolverR.get();

            if (convL != nullptr && convR != nullptr)
           #else
            if (TwoStageThreadedConvolver* const conv = convolver.get())
           #endif
            {
               #ifdef CONVOLUTION_REVERB
                if (buffered)
                {
                    const float* const ins[2] = { highpassBufL, highpassBufR };
                    bufferedConvolver.process(ins, outputs, frames);
                }
                else
               #endif
                {
                   #ifdef CONVOLUTION_REVERB
                    convL->process(highpassBufL, outL, frames);
                    convR->process(highpassBufR, outR, frames);
                   #else
                    conv->process(highpassBufL, outL, frames);
                   #endif
                }

                for (uint32_t i = 0; i < frames; ++i)
                {
                    dryLevel = smoothDryLevel.next();
                    wetLevel = smoothWetLevel.next();

                    if (wetLevel <= 0.001f)
                    {
                       #ifdef CONVOLUTION_REVERB
                        outL[i] = outR[i] = 0.f;
                       #else
                        outL[i] = 0.f;
                       #endif
                    }
                    else
                    {
                        outL[i] *= wetLevel;
                       #ifdef CONVOLUTION_REVERB
                        outR[i] *= wetLevel;
                       #endif
                    }

                    if (dryLevel > 0.001f)
                    {
                        outL[i] += dryBufL[i] * dryLevel;
                       #ifdef CONVOLUTION_REVERB
                        outR[i] += dryBufR[i] * dryLevel;
                       #endif
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
           #ifdef CONVOLUTION_REVERB
            outR[i] = dryBufR[i] * dryLevel;
           #endif
        }
    }

    void bufferSizeChanged(const uint32_t newBufferSize) override
    {
        bufferSize = newBufferSize;

        delete[] highpassBufL;
        delete[] inplaceProcBufL;

        highpassBufL = new float[newBufferSize];
        inplaceProcBufL = new float[newBufferSize];

       #ifdef CONVOLUTION_REVERB
        delete[] highpassBufR;
        delete[] inplaceProcBufR;

        highpassBufR = new float[newBufferSize];
        inplaceProcBufR = new float[newBufferSize];

        bufferedConvolver.setBufferSize(newBufferSize);
       #endif
    }

    void sampleRateChanged(const double newSampleRate) override
    {
       #ifdef CONVOLUTION_REVERB
        korgFilterL.setSampleRate(newSampleRate);
        korgFilterR.setSampleRate(newSampleRate);
       #endif

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
   #ifdef CONVOLUTION_REVERB
    ScopedPointer<TwoStageThreadedConvolver> convolverL;
    ScopedPointer<TwoStageThreadedConvolver> convolverR;
    StereoBufferedConvolver bufferedConvolver;
    Korg35Filter korgFilterL, korgFilterR;
   #else
    ScopedPointer<TwoStageThreadedConvolver> convolver;
   #endif

    Mutex mutex;
    String loadedFilename;

    bool bypassed = false;
   #ifdef CONVOLUTION_REVERB
    bool buffered = false;
    bool trails = true;
   #else
    static constexpr const bool trails = false;
   #endif
    uint32_t bufferSize = 0;

    float parameters[kParameterCount];

    // smoothed parameters
    LinearValueSmoother smoothDryLevel;
    LinearValueSmoother smoothWetLevel;

    // buffers for placing highpass signal before convolution
    float* highpassBufL = nullptr;
   #ifdef CONVOLUTION_REVERB
    float* highpassBufR = nullptr;
   #endif

    // if doing inline processing, copy buffers here before convolution
    float* inplaceProcBufL = nullptr;
   #ifdef CONVOLUTION_REVERB
    float* inplaceProcBufR = nullptr;
   #endif

    DISTRHO_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ConvolutionLoaderPlugin)
};

// -----------------------------------------------------------------------

Plugin* createPlugin()
{
    return new ConvolutionLoaderPlugin();
}

END_NAMESPACE_DISTRHO
