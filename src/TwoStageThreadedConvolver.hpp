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
        if (irLen > headBlockSize * 2)
        {
            if (! fftconvolver::TwoStageFFTConvolver::init(headBlockSize, tailBlockSize, ir, irLen))
                return false;

            startThread(true);
            return true;
        }

        nonThreadedConvolver = new fftconvolver::FFTConvolver();
        return true;
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

END_NAMESPACE_DISTRHO
