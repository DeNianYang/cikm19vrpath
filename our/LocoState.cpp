#include "LocoState.h"

LocoState::LocoState()
{
}

LocoState::LocoState(int v_x, int v_y, int in_theta_v, int p_x, int p_y, int in_theta_p)
{
    gamma_v_x = v_x;
    gamma_v_y = v_y;
    theta_v = in_theta_v % 360;
    gamma_p_x = p_x;
    gamma_p_y = p_y;
    theta_p = in_theta_p % 360;
    valid = true;
}

LocoState::LocoState(const LocoState& lc)
{
    gamma_v_x = lc.gamma_v_x;
    gamma_v_y = lc.gamma_v_y;
    theta_v = lc.theta_v;
    gamma_p_x = lc.gamma_p_x;
    gamma_p_y = lc.gamma_p_y;
    theta_p = lc.theta_p;
    valid = lc.valid;
}

int LocoState::getThetaV() {
    return theta_v;
}

int LocoState::getThetaP() {
    return theta_p;
}

int LocoState::getGammaPX() {
    return gamma_p_x;
}

int LocoState::getGammaPY() {
    return gamma_p_y;
}

int LocoState::getGammaVX() {
    return gamma_v_x;
}

int LocoState::getGammaVY() {
    return gamma_v_y;
}

bool LocoState::getValid() {
    return valid;
}

std::string LocoState::to_string() {
    std::stringstream ss_1;
    ss_1 << gamma_v_x << " " << gamma_v_y << " " << theta_v << " " << gamma_p_x << " " << gamma_p_y << " " << theta_p;
    return ss_1.str();
}

void LocoState::sumPosition(int delta_v_x, int delta_v_y, int delta_p_x, int delta_p_y) {
    gamma_v_x += delta_v_x;
    gamma_v_y += delta_v_y;
    gamma_p_x += delta_p_x;
    gamma_p_y += delta_p_y;
}

void LocoState::sumAngle(int rotate_v, int rotate_p) {
    theta_v = (theta_v + rotate_v + 360) % 360;
    theta_p = (theta_p + rotate_p + 360) % 360;
}

bool LocoState::isValid(int width_v, int length_v, int width_p, int length_p) {
    valid = (0 <= gamma_v_x && gamma_v_x < width_v && 0 <= gamma_v_y && gamma_v_y < length_v
           && 0 <= gamma_p_x && gamma_p_x < width_p && 0 <= gamma_p_y && gamma_p_y < length_p);
    return valid;
}

void LocoState::print_state() {
    printf("gamma_v = (%d, %d) theta_v = %d gamma_p = (%d, %d) theta_p = %d",
        gamma_v_x, gamma_v_y, theta_v, gamma_p_x, gamma_p_y, theta_p);
}

bool LocoState::operator ==(const LocoState &st) const {
    return ((gamma_v_x == st.gamma_v_x) && (gamma_v_y == st.gamma_v_y)
    && (gamma_p_x == st.gamma_p_x) && (gamma_p_y == st.gamma_p_y)
    && (theta_v == st.theta_v) && (theta_p == st.theta_p));        
}

bool LocoState::operator <(const LocoState &st) const {
    
    std::stringstream ss_1;
    ss_1 << gamma_v_x << " " << gamma_v_y << " " << theta_v << " " << gamma_p_x << " " << gamma_p_y << " " << theta_p;
    std::string s_1 = ss_1.str();

    std::stringstream ss_2;
    ss_2 << st.gamma_v_x << " " << st.gamma_v_y << " " << st.theta_v << " " << st.gamma_p_x << " " << st.gamma_p_y << " " << st.theta_p;
    std::string s_2 = ss_2.str();
    return s_1 < s_2;
}

LocoState LocoState::copyState() {
    LocoState ls(gamma_v_x, gamma_v_y, theta_v, gamma_p_x, gamma_p_y, theta_p);
    return ls;
}

void LocoState::print() {
    std::cout << "(" << gamma_v_x << "," << gamma_v_y << "," << theta_v << ", " << gamma_p_x << "," << gamma_p_y << "," << theta_p << ")" << std::endl;
}

LocoState::~LocoState()
{

}
