#ifndef LocoState_H
#define LocoState_H

#include "Tools.h"

class LocoState
{
public:
    LocoState();
    LocoState(int v_x, int v_y, int in_theta_v, int p_x, int p_y, int in_theta_p);
    LocoState(const LocoState& lc);
    ~LocoState();

    int getThetaV();
    int getThetaP();
    int getGammaPX();
    int getGammaPY();
    int getGammaVX();
    int getGammaVY();
    bool getValid();
    void sumPosition(int delta_v_x, int delta_v_y, int delta_p_x, int delta_p_y);
    void sumAngle(int rotate_v, int rotate_p);
    void print_state();
    std::string to_string();
    LocoState copyState();

    
    bool operator ==(const LocoState &st) const;
    bool operator <(const LocoState &st) const;

    void print();
    
protected:

    bool isValid(int width_v, int length_v, int width_p, int length_p);

    int gamma_v_x;
    int gamma_v_y;
    int theta_v;
    int gamma_p_x;
    int gamma_p_y;
    int theta_p;
    bool valid;

};
#endif // LocoState_H
