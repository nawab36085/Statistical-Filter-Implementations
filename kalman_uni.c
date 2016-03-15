#include <stdio.h>
#include <math.h>
float* update(float old_mean, float new_mean, float old_variance, float new_variance);
float* predict(float old_mean, float motion_mean, float old_variance, float motion_variance);
float updated_values[2] = {1, 100000};
float* values = updated_values;
int main(){
        int i;
        float measurement[10] = {10, 19, 29, 38, 51, 59, 69, 80, 91, 101};
        float measurement_variance = 1;
        float motion[10] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
        float motion_variance = 1;
        for(i=0; i<sizeof(measurement)/sizeof(float); i++){
                        values = update(values[0], values[1], measurement[i], measurement_variance);
                        values = predict(values[0], values[1], motion[i], motion_variance );
        }

}


float* update(float old_mean, float old_variance, float new_mean, float new_variance){
        //float values[2];
        updated_values[0] = (old_mean*new_variance + new_mean*old_variance)/(old_variance + new_variance);
        updated_values[1] = 1/(1/old_variance + 1/new_variance);
        printf("updated position=%f and updated variance=%f\n", updated_values[0], updated_values[1]);
        return updated_values;
}

float* predict(float old_mean, float old_variance, float motion_mean, float motion_variance){
        //float predicted_values[2];
        updated_values[0] = old_mean + motion_mean;
        updated_values[1] = old_variance + motion_variance;
        printf("predicted position=%f and predicted variance=%f\n", updated_values[0], updated_values[1]);
        return updated_values;
}

