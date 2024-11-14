#include "simulation_base.hpp"

using namespace cgp;

// Convert a density value to a pressure
static float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness*(rho-rho0);
}

static float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
    //  Laplacian of W_viscosity
    float const r = norm(p_i-p_j);
    assert_cgp_no_msg(r<=h);
    return 45.0f/(3.14159f*std::pow(h,6)) * (h-r);
}

static vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
    //  Gradient of W_spiky
    vec3 const& p = p_i-p_j;
    float const r = norm(p);
    assert_cgp_no_msg(r<=h);
    return -45.0f/(3.14159f*std::pow(h,6)) * (p/r) * std::pow(h-r, 2.0f);
}

static float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i-p_j);
    assert_cgp_no_msg(r<=h);
	return 315.0/(64.0*3.14159f*std::pow(h,9)) * std::pow(h*h-r*r, 3.0f);
}

void update_parameters(spatial_grid_base& grid, sph_parameters_structure const& sph_parameters)
{
    float h = sph_parameters.h;

    // Parcourir chaque cellule de la grille
    #pragma omp parallel for
    for (int i = 0; i < grid.grid.size(); ++i)
    {
        // Récupérer les particules dans la cellule courante
        std::vector<particle_element*>& particles_in_cell = grid.grid[i];

        // Obtenir les particules voisines dans les cellules environnantes (3x3 cellules voisines)
        std::vector<particle_element*> neighbors = grid.get_neighbors(i);

        // Pour chaque particule dans la cellule courante
        for (particle_element* particle : particles_in_cell)
        {
            float m = 0.0;
            float nu = 0.0;
            vec3 color{0.0, 0.0, 0.0};
            float div = 0.0;

            auto fluid_i = particle->fluid_type;
            const auto& p_i = particle->p;
            auto v_i = norm(particle->v);

            // Calculer les paramètres à partir des particules voisines
            for (particle_element* neighbor : neighbors)
            {
                if (neighbor != particle) {
                    auto const& p_j = neighbor->p;
                    float const r = norm(p_i - p_j);

                    if (r < h)
                    {
                        auto fluid_j = neighbor->fluid_type;

                        // Vérifier si les fluides sont compatibles pour le mélange
                        if (fluid_i == fluid_j || fluid_i->soluble_classes.find(fluid_j) != fluid_i->soluble_classes.end())
                        {
                            m += neighbor->m;
                            nu += neighbor->nu;
                            color += neighbor->color;
                            div += 1.0;
                        }
                    }
                }
            }

            // Mise à jour des paramètres si des voisins valides ont été trouvés
            if (div > 0.0)
            {
                float fixed_rate = 0.1f;

                float diff_m = (m / div) - particle->m;
                particle->m += sph_parameters.fluid_mixing_rate * diff_m * fixed_rate * v_i;

                float diff_nu = (nu / div) - particle->nu;
                particle->nu += sph_parameters.fluid_mixing_rate * diff_nu * fixed_rate * v_i;

                vec3 diff_color = (color / div) - particle->color;
                particle->color += sph_parameters.fluid_mixing_rate * diff_color * fixed_rate * v_i;
            }
        }
    }
}

void update_density(spatial_grid_base& grid, float h)
{
    // Parcourir chaque cellule de la grille pour calculer la densité des particules
    #pragma omp parallel for
    for (int i = 0; i < grid.grid.size(); ++i)
    {
        // Récupérer les particules dans la cellule courante
        std::vector<particle_element*>& particles_in_cell = grid.grid[i];

        // Obtenir les particules voisines dans les cellules environnantes (3x3 cellules voisines)
        std::vector<particle_element*> neighbors = grid.get_neighbors(i);

        // Pour chaque particule dans la cellule courante
        for (particle_element* particle : particles_in_cell)
        {
            float rho = 0.0f;
            auto const& p_i = particle->p;

            // Calculer la densité à partir des particules voisines
            for (particle_element* neighbor : neighbors)
            {
                auto const& p_j = neighbor->p;
                float const m_j = neighbor->m;

                float const r = norm(p_i - p_j);
                if (r < h)
                {
                    float const w = W_density(p_i, p_j, h);
                    rho += m_j * w;
                }
            }

            // Mettre à jour la densité de la particule
            particle->rho = rho;
        }
    }
}

// Convert the particle density to pressure
void update_pressure(numarray<particle_element>& particles, float rho0, float stiffness)
{
	const int N = particles.size();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
        particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
}

void update_force(spatial_grid_base& grid, sph_parameters_structure const& sph_parameters)
{
    float h = sph_parameters.h;

    // Parcourir chaque cellule de la grille
    #pragma omp parallel for
    for (int i = 0; i < grid.grid.size(); ++i)
    {
        // Récupérer les particules dans la cellule courante
        std::vector<particle_element*>& particles_in_cell = grid.grid[i];

        // Obtenir les particules voisines dans les cellules environnantes (3x3 cellules voisines)
        std::vector<particle_element*> neighbors = grid.get_neighbors(i);

        // Pour chaque particule dans la cellule courante
        for (particle_element* particle : particles_in_cell)
        {
            float m_i = particle->m;
            particle->f = m_i * vec3{0, -1.0f, 0} * sph_parameters.gravity_strength;

            particle->f += particle->external_forces;
            particle->external_forces = {0, 0, 0};

            vec3 F_pressure{0, 0, 0};
            vec3 F_viscosity{0, 0, 0};

            vec3 const& p_i = particle->p;
            float pr_i = particle->pressure;
            float rho_i = particle->rho;
            vec3 const& v_i = particle->v;
            float nu_i = particle->nu;

            // Calculer les forces à partir des particules voisines
            for (particle_element* neighbor : neighbors)
            {
                if (neighbor != particle) {  // Ne pas comparer la particule avec elle-même
                    vec3 const& p_j = neighbor->p;
                    float pr_j = neighbor->pressure;
                    float rho_j = neighbor->rho;
                    vec3 const& v_j = neighbor->v;
                    float m_j = neighbor->m;
                    float nu_j = neighbor->nu;

                    float const r = norm(p_i - p_j);
                    if (r < h)
                    {
                        // Compute F_pressure
                        vec3 grad = W_gradient_pressure(p_i, p_j, h);
                        F_pressure += m_j * (pr_i + pr_j) / (2 * rho_j) * grad;

                        // Compute F_viscosity
                        float lap = W_laplacian_viscosity(p_i, p_j, h);
                        F_viscosity += ((nu_j + nu_i) / 2) * m_j * (v_j - v_i) / rho_j * lap;
                    }
                }
            }

            // Appliquer les forces de pression et de viscosité
            F_pressure *= -m_i / rho_i;
            F_viscosity *= m_i;

            particle->f += F_pressure + F_viscosity;
        }
    }
}