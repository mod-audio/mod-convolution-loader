/*
 * 2-Stage Threaded Convolver
 * Copyright (C) 2022-2023 Filipe Coelho <falktx@falktx.com>
 * SPDX-License-Identifier: ISC
 */

#pragma once

#include "Semaphore.hpp"
#include "extra/ScopedPointer.hpp"
#include "extra/Thread.hpp"

#include "FFTConvolver/TwoStageFFTConvolver.h"

START_NAMESPACE_DISTRHO

// --------------------------------------------------------------------------------------------------------------------

class TwoStageThreadedConvolver : public fftconvolver::TwoStageFFTConvolver,
                                  private Thread
{
    ScopedPointer<fftconvolver::FFTConvolver> nonThreadedConvolver;
    Semaphore semBgProcStart;
    Semaphore semBgProcFinished;

public:
    TwoStageThreadedConvolver()
        : fftconvolver::TwoStageFFTConvolver(),
          Thread("TwoStageThreadedConvolver"),
          semBgProcStart(1),
          semBgProcFinished(0)
    {
    }

    ~TwoStageThreadedConvolver() override
    {
        if (nonThreadedConvolver != nullptr)
        {
            nonThreadedConvolver = nullptr;
            return;
        }

        signalThreadShouldExit();
        semBgProcStart.post();
        stopThread(5000);
    }

    bool init(const size_t headBlockSize, const size_t tailBlockSize, const fftconvolver::Sample* const ir, const size_t irLen)
    {
        if (fftconvolver::TwoStageFFTConvolver::init(headBlockSize, tailBlockSize, ir, irLen))
        {
            startThread(true);
            return true;
        }

        ScopedPointer<fftconvolver::FFTConvolver> conv(new fftconvolver::FFTConvolver);

        if (conv->init(headBlockSize, ir, irLen))
        {
            nonThreadedConvolver = conv.release();
            return true;
        }

        return false;
    }

    void process(const fftconvolver::Sample* const input, fftconvolver::Sample* const output, const size_t len)
    {
        if (nonThreadedConvolver != nullptr)
            nonThreadedConvolver->process(input, output, len);
        else
            fftconvolver::TwoStageFFTConvolver::process(input, output, len);
    }

protected:
    void startBackgroundProcessing() override
    {
        semBgProcStart.post();
    }

    void waitForBackgroundProcessing() override
    {
        if (isThreadRunning() && !shouldThreadExit())
            semBgProcFinished.wait();
    }

    void run() override
    {
        while (!shouldThreadExit())
        {
            semBgProcStart.wait();

            if (shouldThreadExit())
                break;

            doBackgroundProcessing();
            semBgProcFinished.post();
        }
    }

    DISTRHO_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TwoStageThreadedConvolver)
};

// --------------------------------------------------------------------------------------------------------------------

class MonoBufferedConvolver : private Thread
{
    TwoStageThreadedConvolver* conv = nullptr;
    float* bufferedInput = nullptr;
    float* bufferedOutput = nullptr;
    uint32_t bufferSize = 0;

    Semaphore semBgProcStart;
    Semaphore semBgProcFinished;

    // Mutex mutexI, mutexO;

public:
    MonoBufferedConvolver()
        : Thread("MonoBufferedConvolver"),
          semBgProcStart(0),
          semBgProcFinished(0)
    {
    }

    ~MonoBufferedConvolver() override
    {
        stop();

        delete[] bufferedInput;
        delete[] bufferedOutput;
    }

    void setBufferSize(const uint32_t newBufferSize)
    {
        if (bufferSize == newBufferSize)
            return;

        const bool wasRunning = isThreadRunning();
        TwoStageThreadedConvolver* const c = conv;

        if (wasRunning)
            stop();

        bufferSize = newBufferSize;

        delete[] bufferedInput;
        delete[] bufferedOutput;

        bufferedInput = new float[newBufferSize];
        bufferedOutput = new float[newBufferSize];
        std::memset(bufferedOutput, 0, sizeof(float)*newBufferSize);

        if (wasRunning)
            start(c);
    }

    void start(TwoStageThreadedConvolver* const c)
    {
        DISTRHO_SAFE_ASSERT_RETURN(c != nullptr,);

        conv = c;
        startThread(true);
    }

    void stop()
    {
        if (conv == nullptr)
            return;

        conv = nullptr;
        signalThreadShouldExit();
        semBgProcStart.post();
        stopThread(5000);
    }

    void process(const float* const input, float* const output, const uint32_t len)
    {
        if (conv == nullptr)
        {
            if (output != input)
                std::memcpy(output, input, sizeof(float)*len);
            return;
        }

        DISTRHO_SAFE_ASSERT_UINT2_RETURN(bufferSize == len, bufferSize, len,);

        // fetch audio from previous processing
        semBgProcFinished.wait();
        {
            // const MutexLocker cmtl(mutexO);
            std::memcpy(output, bufferedOutput, sizeof(float)*len);
        }

        // place input for next processing
        {
            // const MutexLocker cmtl(mutexI);
            std::memcpy(bufferedInput, input, sizeof(float)*len);
        }
        semBgProcStart.post();
    }

    void run() override
    {
        DISTRHO_SAFE_ASSERT_RETURN(bufferSize != 0,);

        TwoStageThreadedConvolver* c;

        float* tmpInput = new float[bufferSize];
        float* tmpOutput = new float[bufferSize];

        std::memset(tmpOutput, 0, sizeof(float)*bufferSize);

        while (!shouldThreadExit())
        {
            semBgProcFinished.post();
            semBgProcStart.wait();

            if (shouldThreadExit())
                break;

            c = conv;
            DISTRHO_SAFE_ASSERT_CONTINUE(c != nullptr);

            {
                // const MutexLocker cmtl(mutexI);
                std::memcpy(tmpInput, bufferedInput, sizeof(float)*bufferSize);
            }

            c->process(tmpInput, tmpOutput, bufferSize);

            {
                // const MutexLocker cmtl(mutexO);
                std::memcpy(bufferedOutput, tmpOutput, sizeof(float)*bufferSize);
            }
        }

        delete[] tmpInput;
        delete[] tmpOutput;
    }

    DISTRHO_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MonoBufferedConvolver)
};

// --------------------------------------------------------------------------------------------------------------------

class StereoBufferedConvolver : private Thread
{
    TwoStageThreadedConvolver* convL = nullptr;
    TwoStageThreadedConvolver* convR = nullptr;
    float* bufferedInputL = nullptr;
    float* bufferedInputR = nullptr;
    float* bufferedOutputL = nullptr;
    float* bufferedOutputR = nullptr;
    uint32_t bufferSize = 0;

    Semaphore semBgProcStart;
    Semaphore semBgProcFinished;

    // Mutex mutexI, mutexO;

public:
    StereoBufferedConvolver()
        : Thread("StereoBufferedConvolver"),
          semBgProcStart(0),
          semBgProcFinished(0)
    {
    }

    ~StereoBufferedConvolver() override
    {
        stop();

        delete[] bufferedInputL;
        delete[] bufferedInputR;
        delete[] bufferedOutputL;
        delete[] bufferedOutputR;
    }

    void setBufferSize(const uint32_t newBufferSize)
    {
        if (bufferSize == newBufferSize)
            return;

        const bool wasRunning = isThreadRunning();
        TwoStageThreadedConvolver* const cl = convL;
        TwoStageThreadedConvolver* const cr = convR;

        if (wasRunning)
            stop();

        bufferSize = newBufferSize;

        delete[] bufferedInputL;
        delete[] bufferedInputR;
        delete[] bufferedOutputL;
        delete[] bufferedOutputR;

        bufferedInputL = new float[newBufferSize];
        bufferedInputR = new float[newBufferSize];
        bufferedOutputL = new float[newBufferSize];
        bufferedOutputR = new float[newBufferSize];
        std::memset(bufferedOutputL, 0, sizeof(float)*newBufferSize);
        std::memset(bufferedOutputR, 0, sizeof(float)*newBufferSize);

        if (wasRunning)
            start(cl, cr);
    }

    void start(TwoStageThreadedConvolver* const cl, TwoStageThreadedConvolver* const cr)
    {
        DISTRHO_SAFE_ASSERT_RETURN(cl != nullptr,);
        DISTRHO_SAFE_ASSERT_RETURN(cr != nullptr,);

        convL = cl;
        convR = cr;
        startThread(true);
    }

    void stop()
    {
        if (convL == nullptr)
            return;

        convL = nullptr;
        convR = nullptr;
        signalThreadShouldExit();
        semBgProcStart.post();
        stopThread(5000);
    }

    void process(const float* const inputs[2], float* const outputs[2], const uint32_t len)
    {
        if (convL == nullptr || convR == nullptr)
        {
            if (outputs[0] != inputs[0])
                std::memcpy(outputs[0], inputs[0], sizeof(float)*len);
            if (outputs[1] != inputs[1])
                std::memcpy(outputs[1], inputs[1], sizeof(float)*len);
            return;
        }

        DISTRHO_SAFE_ASSERT_UINT2_RETURN(bufferSize == len, bufferSize, len,);

        // fetch audio from previous processing
        semBgProcFinished.wait();
        {
            // const MutexLocker cmtl(mutexO);
            std::memcpy(outputs[0], bufferedOutputL, sizeof(float)*len);
            std::memcpy(outputs[1], bufferedOutputR, sizeof(float)*len);
        }

        // place input for next processing
        {
            // const MutexLocker cmtl(mutexI);
            std::memcpy(bufferedInputL, inputs[0], sizeof(float)*len);
            std::memcpy(bufferedInputR, inputs[1], sizeof(float)*len);
        }
        semBgProcStart.post();
    }

    void run() override
    {
        DISTRHO_SAFE_ASSERT_RETURN(bufferSize != 0,);

        TwoStageThreadedConvolver* cl;
        TwoStageThreadedConvolver* cr;

        float* tmpInputL = new float[bufferSize];
        float* tmpInputR = new float[bufferSize];
        float* tmpOutputL = new float[bufferSize];
        float* tmpOutputR = new float[bufferSize];

        std::memset(tmpOutputL, 0, sizeof(float)*bufferSize);
        std::memset(tmpOutputR, 0, sizeof(float)*bufferSize);

        while (!shouldThreadExit())
        {
            semBgProcFinished.post();
            semBgProcStart.wait();

            if (shouldThreadExit())
                break;

            cl = convL;
            cr = convR;
            DISTRHO_SAFE_ASSERT_CONTINUE(cl != nullptr);
            DISTRHO_SAFE_ASSERT_CONTINUE(cr != nullptr);

            {
                // const MutexLocker cmtl(mutexI);
                std::memcpy(tmpInputL, bufferedInputL, sizeof(float)*bufferSize);
                std::memcpy(tmpInputR, bufferedInputR, sizeof(float)*bufferSize);
            }

            cl->process(tmpInputL, tmpOutputL, bufferSize);
            cr->process(tmpInputR, tmpOutputR, bufferSize);

            {
                // const MutexLocker cmtl(mutexO);
                std::memcpy(bufferedOutputL, tmpOutputL, sizeof(float)*bufferSize);
                std::memcpy(bufferedOutputR, tmpOutputR, sizeof(float)*bufferSize);
            }
        }

        delete[] tmpInputL;
        delete[] tmpInputR;
        delete[] tmpOutputL;
        delete[] tmpOutputR;
    }

    DISTRHO_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(StereoBufferedConvolver)
};

// --------------------------------------------------------------------------------------------------------------------

END_NAMESPACE_DISTRHO
